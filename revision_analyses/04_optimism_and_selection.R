# ==============================================================================
# 04. Bootstrap optimism correction and stability of variable selection
#
# Reproduces (manuscript numbering):
#   Supplementary Table 6  optimism-corrected performance, two versions
#   Supplementary Table 7  stability of predictor selection (full-pipeline version)
#
# 100 bootstrap replicates of the derivation partition (Harrell's optimism correction), in two versions:
#   "variables fixed"   the seven published predictors are held fixed; cutoffs, weights and the
#                       score-to-probability map are re-derived in every replicate;
#   "full pipeline"     variable selection is repeated too, by an AUTOMATED SIMPLIFICATION of the selection rule:
#                       the seven highest-ranked variables by XGBoost information gain among the eligible candidates.
#                       The redundancy criterion and clinical interpretability of the actual rule cannot be resampled.
# Eligible candidates: numeric or binary variables obtainable at triage (chief complaints, past medical history,
# laboratory, imaging, previous disposition, arrival and demographic/administrative variables are excluded).
# Missing values are imputed with the median of each bootstrap sample (repeating multiple imputation 100 times is
# not feasible), which makes the correction slightly anti-conservative with respect to imputation uncertainty.
#
# The replicates are independent, so they run in parallel (set PANACEA_CORES); results do not depend on the number
# of cores because the bootstrap indices are drawn first, in a fixed order.
# From the repository root:  Rscript revision_analyses/04_optimism_and_selection.R   (about 6 minutes on 8 cores)
# ==============================================================================
source("revision_analyses/R/setup.R")
library(parallel)
B <- 100; RANK_NROUNDS <- 50; N_ITEMS <- 7
N_CORES <- as.integer(Sys.getenv("PANACEA_CORES", max(1, detectCores() - 2)))
PUBLISHED_ORDER <- c("esi", "age", "meds_cardiovascular", "meds_analgesics", "meds_psychotherapeuticdrugs", "n_admissions", "triage_vital_hr")

# ---- eligible candidates ---------------------------------------------------------
exclude_prefix <- c("^cc_", "^pmh_", "^lab", "^img", "^previousdispo", "^arrival", "^dep_name", "^gender", "^race", "^ethnic",
                    "^lang", "^religion", "^marital", "^employ", "^insurance")
pool <- setdiff(names(train_data_raw), c("id", "disposition", "disposition_bin"))
for (px in exclude_prefix) pool <- pool[!grepl(px, pool)]
as_num <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x)
usable <- sapply(pool, function(v) { x <- suppressWarnings(as_num(train_data_raw[[v]])); !all(is.na(x)) && length(unique(na.omit(x))) >= 2 })
POOL_VARS <- union(pool[usable], PUBLISHED_ORDER)
cat("Eligible candidates:", length(POOL_VARS), "\n")
work <- as.data.frame(lapply(train_data_raw[, POOL_VARS, drop = FALSE], as_num)); names(work) <- POOL_VARS
work$disposition_bin <- y_train
impute_within <- function(d) { for (v in setdiff(names(d), "disposition_bin")) if (anyNA(d[[v]])) d[[v]][is.na(d[[v]])] <- median(d[[v]], na.rm = TRUE); d }
work <- impute_within(work)

# ---- pipeline components -----------------------------------------------------------
rank_variables <- function(data) {
  vars <- setdiff(names(data), "disposition_bin"); X <- as.matrix(data[, vars, drop = FALSE]); storage.mode(X) <- "double"
  fit <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 6, eta = 0.3, nthread = 2),
                   data = xgb.DMatrix(X, label = data$disposition_bin), nrounds = RANK_NROUNDS, verbose = 0)
  imp <- xgb.importance(model = fit); imp$Feature[order(-imp$Gain)]
}
youden_cut <- function(y, x) as.numeric(pROC::coords(pROC::roc(y, x, quiet = TRUE), "best", best.method = "youden", ret = "threshold")[1])
build_score <- function(data, vars, cutoffs, weights) {
  s <- rep(0, nrow(data))
  for (v in vars) { w <- weights[[v]]
    if (v == "esi") s <- s + w * (5 - data$esi) else if (v == "age") s <- s + w * pmin(floor(data$age / 10), 9) else s <- s + ifelse(data[[v]] > cutoffs[[v]], w, 0) }
  s
}
optimize_weights <- function(data, vars, cutoffs) {                      # greedy coordinate ascent, explicit order, weights 1-5, one pass
  w <- setNames(as.list(rep(1, length(vars))), vars); y <- data$disposition_bin
  for (v in vars) { best_w <- w[[v]]; best_a <- -Inf
    for (cand in 1:5) { w[[v]] <- cand; a <- fast_auc(y, build_score(data, vars, cutoffs, w)); if (!is.na(a) && a > best_a) { best_a <- a; best_w <- cand } }
    w[[v]] <- best_w }
  w
}
derive_pipeline <- function(data, mode) {
  vars <- if (mode == "variables_fixed") PUBLISHED_ORDER else head(rank_variables(data), N_ITEMS)
  cont <- setdiff(vars, c("esi", "age"))
  cutoffs <- setNames(lapply(cont, function(v) youden_cut(data$disposition_bin, data[[v]])), cont)
  w <- optimize_weights(data, vars, cutoffs); s <- build_score(data, vars, cutoffs, w)
  mp <- data.frame(s = s, y = data$disposition_bin) %>% mutate(g = pmin(pmax(s, 4), 33)) %>% group_by(g) %>% summarise(p = mean(y), .groups = "drop")
  list(vars = vars, cutoffs = cutoffs, weights = w, map = mp)
}
apply_model <- function(m, data) {
  g <- pmin(pmax(build_score(data, m$vars, m$cutoffs, m$weights), 4), 33); p <- m$map$p[match(g, m$map$g)]
  if (anyNA(p)) { i <- which(is.na(p)); p[i] <- m$map$p[sapply(g[i], function(z) which.min(abs(m$map$g - z)))] }
  p
}
evaluate <- function(p, y) {
  a <- fast_auc(y, p); pc <- pmin(pmax(p, 1e-6), 1 - 1e-6); lp <- qlogis(pc)
  c(AUC = a, Slope = unname(coef(glm(y ~ lp, family = binomial))[2]), Intercept = unname(coef(glm(y ~ offset(lp), family = binomial))[1]),
    Brier = mean((pc - y)^2), OE = sum(y) / sum(pc))
}

# ---- apparent performance of the published pipeline -----------------------------------
m_pub <- derive_pipeline(work, "variables_fixed")
stopifnot(all(unlist(m_pub$weights) == PUBLISHED_WEIGHTS[PUBLISHED_ORDER]))
perf_pub <- evaluate(apply_model(m_pub, work), work$disposition_bin)
cat("Apparent performance:\n"); print(round(perf_pub, 5))

# ---- Harrell's optimism correction -------------------------------------------------------
run_optimism <- function(mode) {
  set.seed(20260813)                                                    # same bootstrap draws for both versions
  idx_list <- lapply(seq_len(B), function(b) sample.int(nrow(work), nrow(work), replace = TRUE))
  res <- mclapply(seq_len(B), function(b) {
    boot <- impute_within(work[idx_list[[b]], ])
    m <- tryCatch(derive_pipeline(boot, mode), error = function(e) NULL); if (is.null(m)) return(NULL)
    list(opt = evaluate(apply_model(m, boot), boot$disposition_bin) - evaluate(apply_model(m, work), work$disposition_bin), vars = m$vars, w = unlist(m$weights))
  }, mc.cores = N_CORES)
  res <- Filter(Negate(is.null), res)
  list(optimism = colMeans(do.call(rbind, lapply(res, `[[`, "opt"))), selected = unlist(lapply(res, `[[`, "vars")), weights = do.call(rbind, lapply(res, `[[`, "w")))
}
res_fixed <- run_optimism("variables_fixed"); res_full <- run_optimism("full_pipeline")

corrected <- rbind(data.frame(Analysis = "Apparent", t(perf_pub)),
                   data.frame(Analysis = "Optimism (full pipeline)", t(res_full$optimism)), data.frame(Analysis = "Corrected (full pipeline)", t(perf_pub - res_full$optimism)),
                   data.frame(Analysis = "Optimism (variables fixed)", t(res_fixed$optimism)), data.frame(Analysis = "Corrected (variables fixed)", t(perf_pub - res_fixed$optimism)))
write.csv(corrected, out("supp_table6_optimism_corrected.csv"), row.names = FALSE); print(corrected, digits = 5, row.names = FALSE)
stab <- as.data.frame(table(res_full$selected)) %>% setNames(c("Variable", "Times_selected")) %>% mutate(Pct = 100 * Times_selected / B) %>% arrange(desc(Times_selected))
write.csv(stab, out("supp_table7_selection_stability.csv"), row.names = FALSE); print(stab, row.names = FALSE)
write.csv(as.data.frame(res_fixed$weights), out("weight_stability_variables_fixed.csv"), row.names = FALSE)
cat("Done: outputs in", OUT_DIR, "\n")
