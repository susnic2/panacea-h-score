# ==============================================================================
# 05. Participant flow, missing data and sensitivity analyses
#
# Reproduces (manuscript numbering):
#   Supplementary Table 1   participant flow
#   Supplementary Table 4   sensitivity to alternative predictor sets (gastrointestinal medications added;
#                           disposition of the previous visit in place of prior-year admissions)
#   Supplementary Table 5   mapping of the 48 medication variables, with the three score items marked
#   Supplementary Table 8   missing data by variable
#   Supplementary Table 13  gain, permutation importance and SHAP        Supplementary Figure 1 (ranks)
#   Supplementary Table 14  sensitivity to the heart-rate cutoff and to the encoding of age
#   Supplementary Table 26  medication missingness by arm and the triage-only score
#   Methods                 variable ranking with and without imputation before XGBoost training
#
# From the repository root:  Rscript revision_analyses/05_sensitivity_analyses.R   (about 10 minutes)
# ==============================================================================
source("revision_analyses/R/setup.R")
library(dplyr)

# ---- A. Participant flow (Supplementary Table 1) -------------------------------------
flow <- data.frame(Step = c("Visits in source dataset", "Excluded: disposition missing", "Analysed cohort", "Derivation partition (70%)", "  of which admitted",
                            "  of which discharged", "Validation partition (30%)", "  of which admitted", "  of which discharged"),
                   N = c(N_SOURCE, N_DISPOSITION_MISSING, nrow(df), nrow(train_data_raw), sum(train_data_raw$disposition == "Admit"), sum(train_data_raw$disposition == "Discharge"),
                         nrow(test_data_raw), sum(test_data_raw$disposition == "Admit"), sum(test_data_raw$disposition == "Discharge")))
flow$Pct <- round(100 * flow$N / nrow(df), 2)
write.csv(flow, out("supp_table1_participant_flow.csv"), row.names = FALSE)

# ---- B. Missing data by variable (Supplementary Table 8) --------------------------------
t1_vars <- intersect(c("age", "gender", "ethnicity", "race", "lang", "religion", "maritalstatus", "employstatus", "insurance_status", "dep_name", "arrivalmode",
                       "arrivalmonth", "arrivalday", "arrivalhour_bin", "triage_vital_hr", "triage_vital_sbp", "triage_vital_dbp", "triage_vital_rr", "triage_vital_o2",
                       "triage_vital_o2_device", "triage_vital_temp", "esi", "previousdispo", "n_edvisits", "n_admissions", "n_surgeries", "meds_cardiovascular",
                       "meds_analgesics", "meds_psychotherapeuticdrugs", "meds_gastrointestinal"), names(df))
adm <- df$disposition == "Admit"
miss_tab <- bind_rows(lapply(t1_vars, function(v) { x <- df[[v]]
  data.frame(Variable = v, N_total = length(x), N_missing = sum(is.na(x)), Pct_missing = 100 * mean(is.na(x)),
             Pct_missing_admitted = 100 * mean(is.na(x[adm])), Pct_missing_discharged = 100 * mean(is.na(x[!adm]))) }))
write.csv(miss_tab, out("supp_table8_missingness_by_variable.csv"), row.names = FALSE)

# ---- C. Medication variables (Supplementary Table 5) -------------------------------------
meds_all <- grep("^meds_", names(df), value = TRUE); stopifnot(length(meds_all) == 48)
narrow <- c("meds_cardiovascular", "meds_analgesics", "meds_psychotherapeuticdrugs")             # the three variables that feed the score
meds_tab <- data.frame(Variable = meds_all, Pct_nonzero = sapply(meds_all, function(v) 100 * mean(as.numeric(as.character(df[[v]])) > 0, na.rm = TRUE)), In_score = meds_all %in% narrow)
write.csv(meds_tab, out("supp_table5_medication_variables.csv"), row.names = FALSE)

# ---- D. Alternative predictor sets and the triage-only score (Supplementary Tables 4 and 26) ------
# Each variant repeats the whole pipeline on its own variables: multiple imputation on the derivation partition,
# the cutoffs of the main score, greedy weights 1-5, median/mode imputation of the held-out partition.
score_flex <- function(data, weights, cutoffs) {
  s <- numeric(nrow(data))
  if ("esi" %in% names(weights)) s <- s + (5 - as.numeric(as.character(data$esi))) * weights[["esi"]]
  if ("age" %in% names(weights)) s <- s + pmin(floor(data$age / 10), 9) * weights[["age"]]
  for (med in names(weights)[grepl("^meds_", names(weights))]) s <- s + ifelse(data[[med]] > 0, weights[[med]], 0)
  if ("n_admissions" %in% names(weights)) s <- s + ifelse(data$n_admissions > cutoffs[["n_admissions"]], weights[["n_admissions"]], 0)
  if ("triage_vital_hr" %in% names(weights)) s <- s + ifelse(data$triage_vital_hr > cutoffs[["triage_vital_hr"]], weights[["triage_vital_hr"]], 0)
  if ("previousdispo" %in% names(weights)) {
    lev <- levels(factor(data$previousdispo)); pos <- lev[grepl("^Admi", lev)][1]; stopifnot(!is.na(pos))   # level meaning "previous visit ended in admission"
    s <- s + ifelse(as.character(data$previousdispo) == pos, weights[["previousdispo"]], 0)
  }
  s
}
run_variant <- function(vars) {
  tr <- train_data_raw[, c(vars, "disposition", "id")]
  fit <- mice(tr[, vars], m = 5, maxit = 10, method = "pmm", seed = SEED_MICE, printFlag = FALSE)
  tri <- complete(fit, action = 1); tri$disposition_bin <- y_train
  w <- setNames(as.list(rep(1, length(vars))), vars)
  for (v in vars) { best_w <- 1; best_a <- -Inf
    for (k in 1:5) { cand <- w; cand[[v]] <- k; s <- score_flex(tri, cand, cutoffs); a <- if (length(unique(s)) > 1) fast_auc(y_train, s) else 0.5; if (a > best_a) { best_a <- a; best_w <- k } }
    w[[v]] <- best_w }
  te <- test_data_raw[, c(vars, "disposition", "id")]
  for (v in vars) if (anyNA(te[[v]])) te[[v]][is.na(te[[v]])] <- if (is.numeric(tri[[v]])) median(tri[[v]], na.rm = TRUE) else names(sort(table(tri[[v]]), decreasing = TRUE))[1]
  s <- score_flex(te, w, cutoffs)
  list(auc = fast_auc(y_test, s), delong_p = roc.test(roc(y_test, test_final$score, quiet = TRUE), roc(y_test, s, quiet = TRUE), method = "delong")$p.value, weights = unlist(w))
}
v_gi <- run_variant(c(score_variables, "meds_gastrointestinal"))
v_prev <- run_variant(c(setdiff(score_variables, "n_admissions"), "previousdispo"))
v_tri <- run_variant(c("esi", "age", "n_admissions", "triage_vital_hr"))
auc_main <- fast_auc(y_test, test_final$score)
write.csv(data.frame(Variant = c("Published score (7 items)", "Published score plus gastrointestinal medications (8 items)",
                                 "Published score with disposition of the previous visit in place of prior-year admissions (7 items)"),
                     AUC_test = c(auc_main, v_gi$auc, v_prev$auc), DeLong_p_vs_published_score = c(NA, v_gi$delong_p, v_prev$delong_p)),
          out("supp_table4_predictor_set_sensitivity.csv"), row.names = FALSE)
miss_arm <- train_data_raw %>% group_by(disposition) %>%
  summarise(n = n(), across(all_of(c("meds_cardiovascular", "meds_gastrointestinal", "meds_analgesics", "meds_psychotherapeuticdrugs")),
                            list(pct_missing = ~ round(mean(is.na(.)) * 100, 2), pct_nonzero_amongNonMissing = ~ round(mean(. > 0, na.rm = TRUE) * 100, 2)), .names = "{.col}__{.fn}"))
write.csv(miss_arm, out("supp_table26_panelA_medication_missingness_by_arm.csv"), row.names = FALSE)
write.csv(data.frame(Variant = c("Published score (7 items)", "Triage-only score (4 items, no medications)"), AUC_test = c(auc_main, v_tri$auc), DeLong_p = c(NA, v_tri$delong_p)),
          out("supp_table26_panelB_triage_only_score.csv"), row.names = FALSE)

# ---- E. Heart-rate cutoff and encoding of age (Supplementary Table 14) -----------------------------
# Heart rate: everything else held at its published weight and cutoff; derivation partition imputed with its median/mode.
imp_sc <- impute_median_mode(train_data_raw[, c(score_variables, "disposition", "id")], test_data_raw[, c(score_variables, "disposition", "id")], vars = score_variables)
sc_train <- imp_sc$train; sc_test <- imp_sc$test
pub_w <- as.list(PUBLISHED_WEIGHTS); pub_cut <- list(n_admissions = 0.5, triage_vital_hr = 98.83)
score_hr <- function(d, w, h) { cc <- pub_cut; cc[["triage_vital_hr"]] <- h; calculate_score(d, w, cc) }
hr_grid <- sort(unique(c(90, 95, 98.83, 100, 105, 110, quantile(sc_train$triage_vital_hr, c(.5, .6, .7, .75, .8, .85, .9)))))
hr_sens <- bind_rows(lapply(hr_grid, function(h) data.frame(HR_cutoff = h, AUC_train = fast_auc(y_train, score_hr(sc_train, pub_w, h)), AUC_test = fast_auc(y_test, score_hr(sc_test, pub_w, h)),
                                                             Pct_train_above_cutoff = mean(sc_train$triage_vital_hr > h))))
hr_sens$Delta_AUC_test_vs_published <- hr_sens$AUC_test - hr_sens$AUC_test[hr_sens$HR_cutoff == 98.83]
reopt <- bind_rows(lapply(c(90, 95, 98.83, 100, 105, 110), function(h) {
  a <- sapply(1:5, function(k) { w <- pub_w; w[["triage_vital_hr"]] <- k; fast_auc(y_train, score_hr(sc_train, w, h)) }); data.frame(cutoff = h, weight = which.max(a), auc_train = max(a)) }))
write.csv(hr_sens, out("supp_table14_hr_cutoff_sensitivity.csv"), row.names = FALSE); write.csv(reopt, out("supp_table14_hr_weight_reoptimised.csv"), row.names = FALSE)
# Age: both encodings re-derived in the same pipeline as the published score (multiple-imputed derivation partition, greedy weights 1-5)
age_cut <- find_youden_cutoff(train_imputed$age, train_imputed$disposition_bin)
score_age_cut <- function(d, w) { s <- (5 - as.numeric(as.character(d$esi))) * w[["esi"]] + ifelse(d$age > age_cut, w[["age"]], 0)
  for (m in c("meds_cardiovascular", "meds_analgesics", "meds_psychotherapeuticdrugs")) s <- s + ifelse(d[[m]] > 0, w[[m]], 0)
  s + ifelse(d$n_admissions > cutoffs[["n_admissions"]], w[["n_admissions"]], 0) + ifelse(d$triage_vital_hr > cutoffs[["triage_vital_hr"]], w[["triage_vital_hr"]], 0) }
wB <- setNames(as.list(rep(1, length(score_variables))), score_variables)
for (v in score_variables) { best_w <- 1; best_a <- -Inf
  for (k in 1:5) { cand <- wB; cand[[v]] <- k; a <- fast_auc(y_train, score_age_cut(train_imputed, cand)); if (a > best_a) { best_a <- a; best_w <- k } }
  wB[[v]] <- best_w }
age_tab <- data.frame(Encoding = c("Decade bands (published)", "Single Youden-optimal cutoff"), Age_cutoff = c(NA, age_cut),
                      Test_AUC = c(auc_main, fast_auc(y_test, score_age_cut(test_final, wB))))
age_tab$Difference <- age_tab$Test_AUC - age_tab$Test_AUC[1]
write.csv(age_tab, out("supp_table14_age_encoding.csv"), row.names = FALSE)

# ---- F. Gain, permutation importance and SHAP (Supplementary Table 13, Supplementary Figure 1) -------------
# Same eligible pool and XGBoost architecture as the automated selection of script 04 (depth 6, eta 0.3, 100 rounds).
exclude_prefix <- c("^cc_", "^pmh_", "^lab", "^img", "^previousdispo", "^arrival", "^dep_name", "^gender", "^race", "^ethnic", "^lang", "^religion", "^marital", "^employ", "^insurance")
pool <- setdiff(names(train_data_raw), c("id", "disposition")); for (px in exclude_prefix) pool <- pool[!grepl(px, pool)]
as_num <- function(x) if (is.factor(x)) as.numeric(as.character(x)) else if (is.character(x)) suppressWarnings(as.numeric(x)) else as.numeric(x)
POOL_VARS <- union(pool[sapply(pool, function(v) !all(is.na(suppressWarnings(as_num(train_data_raw[[v]])))))], score_variables)
cat("Eligible candidates for importance:", length(POOL_VARS), "\n")
wtr <- as.data.frame(lapply(train_data_raw[, POOL_VARS, drop = FALSE], as_num)); names(wtr) <- POOL_VARS
wte <- as.data.frame(lapply(test_data_raw[, POOL_VARS, drop = FALSE], as_num)); names(wte) <- POOL_VARS
impute_within <- function(d, ref = d) { for (v in names(d)) if (anyNA(d[[v]])) d[[v]][is.na(d[[v]])] <- median(ref[[v]], na.rm = TRUE); d }
Xtr <- as.matrix(impute_within(wtr)); Xte <- as.matrix(impute_within(wte, wtr)); storage.mode(Xtr) <- "double"; storage.mode(Xte) <- "double"
set.seed(2024)
xgb_pool <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 6, eta = 0.3, nthread = 2), data = xgb.DMatrix(Xtr, label = y_train), nrounds = 100, verbose = 0)
dte <- xgb.DMatrix(Xte, label = y_test); base_auc <- fast_auc(y_test, predict(xgb_pool, dte))
imp_gain <- xgb.importance(model = xgb_pool)
set.seed(20260813); N_PERM <- 10; PERM_TOPN <- 40
compare_vars <- intersect(c(score_variables, "meds_gastrointestinal"), POOL_VARS)
perm_vars <- union(head(imp_gain$Feature, PERM_TOPN), compare_vars)
perm_imp <- sapply(perm_vars, function(v) { ci <- match(v, colnames(Xte)); mean(sapply(seq_len(N_PERM), function(r) { Xp <- Xte; Xp[, ci] <- sample(Xp[, ci]); base_auc - fast_auc(y_test, predict(xgb_pool, xgb.DMatrix(Xp))) })) })
perm_tab <- data.frame(Feature = perm_vars, Perm_AUC_drop = as.numeric(perm_imp)) %>% arrange(desc(Perm_AUC_drop))
shap <- predict(xgb_pool, dte, predcontrib = TRUE); shap <- shap[, -ncol(shap), drop = FALSE]   # the last column is the intercept ("BIAS" or "(Intercept)", depending on the xgboost version)
shap_tab <- data.frame(Feature = colnames(shap), Mean_abs_SHAP = colMeans(abs(shap))) %>% arrange(desc(Mean_abs_SHAP))
rg <- setNames(seq_len(nrow(imp_gain)), imp_gain$Feature); rp <- setNames(seq_len(nrow(perm_tab)), perm_tab$Feature); rs <- setNames(seq_len(nrow(shap_tab)), shap_tab$Feature)
comparison <- data.frame(Variable = compare_vars, Rank_gain = rg[compare_vars], Rank_permutation = rp[compare_vars], Rank_SHAP = rs[compare_vars],
                         Gain = imp_gain$Gain[match(compare_vars, imp_gain$Feature)], Perm_AUC_drop = perm_tab$Perm_AUC_drop[match(compare_vars, perm_tab$Feature)],
                         Mean_abs_SHAP = shap_tab$Mean_abs_SHAP[match(compare_vars, shap_tab$Feature)], In_published_score = compare_vars %in% score_variables)
write.csv(comparison, out("supp_table13_importance_comparison.csv"), row.names = FALSE)
rk <- data.frame(gain = rg[perm_vars], perm = rp[perm_vars], shap = rs[perm_vars])
cat(sprintf("Spearman: gain-permutation %.3f | gain-SHAP %.3f | permutation-SHAP %.3f\n", cor(rk$gain, rk$perm, method = "spearman"), cor(rk$gain, rk$shap, method = "spearman"), cor(rk$perm, rk$shap, method = "spearman")))

# ---- G. Variable ranking with and without imputation before XGBoost training (Methods) -----------------
set.seed(456); xgb_folds <- createFolds(train_data_raw$disposition, k = 5, list = TRUE, returnTrain = FALSE)
run_ranking <- function(impute) {
  imps <- list(); aucs <- numeric(5)
  for (i in 1:5) { te <- xgb_folds[[i]]; ftr <- train_data_raw[-te, ]; fte <- train_data_raw[te, ]
    if (impute) { im <- impute_median_mode(ftr, fte); ftr <- im$train; fte <- im$test }
    a <- prepare_xgb_data(ftr); b <- prepare_xgb_data(fte)
    xa <- as.matrix(a[, setdiff(names(a), "disposition")]); xb <- as.matrix(b[, setdiff(names(b), "disposition")])
    m <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1), data = xgb.DMatrix(xa, label = a$disposition), nrounds = 100, verbose = 0)
    aucs[i] <- fast_auc(b$disposition, predict(m, xgb.DMatrix(xb, label = b$disposition)))
    g <- xgb.importance(feature_names = colnames(xa), model = m); g$Fold <- i; imps[[i]] <- g }
  list(auc = aucs, rank = bind_rows(imps) %>% group_by(Feature) %>% summarise(Mean_Gain = mean(Gain), .groups = "drop") %>% arrange(desc(Mean_Gain)))
}
A <- run_ranking(TRUE); B <- run_ranking(FALSE)
cmp <- merge(transform(A$rank, Rank_imputed = seq_len(nrow(A$rank)), Gain_imputed = Mean_Gain)[, c("Feature", "Rank_imputed", "Gain_imputed")],
             transform(B$rank, Rank_native = seq_len(nrow(B$rank)), Gain_native = Mean_Gain)[, c("Feature", "Rank_native", "Gain_native")], by = "Feature", all = TRUE)
write.csv(cmp[order(cmp$Rank_imputed), ], out("ranking_with_and_without_imputation.csv"), row.names = FALSE)
write.csv(data.frame(Fold = c(paste0("Fold_", 1:5), "Pooled_Mean"), AUC_imputed = c(A$auc, mean(A$auc)), AUC_native = c(B$auc, mean(B$auc))), out("ranking_cv_auc_with_and_without_imputation.csv"), row.names = FALSE)
cat(sprintf("Pooled 5-fold AUC: imputed %.4f | native NA %.4f | top-15 overlap %d of 15\n", mean(A$auc), mean(B$auc), length(intersect(head(A$rank$Feature, 15), head(B$rank$Feature, 15)))))
cat("Done: outputs in", OUT_DIR, "\n")
