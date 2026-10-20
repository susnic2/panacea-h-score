# ==============================================================================
# Shared setup for the analyses added during peer review
#
# Run every script of this folder from the REPOSITORY ROOT, e.g.
#   Rscript revision_analyses/01_mapping_calibration_bands.R
#
# This file
#   1. loads and assembles the dataset and makes the 70/30 split exactly as
#      Sections 1-3 of panacea_h_reproducibility.R do (it duplicates that code on
#      purpose, so that the original pipeline stays untouched);
#   2. re-derives the PAnACEA-H score on the derivation partition (multiple
#      imputation, Youden cutoffs, greedy weights) and checks that the weights
#      are the published ones;
#   3. builds the held-out partition with the score and caches the objects that
#      take long to compute.
# ==============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(caret); library(xgboost)
  library(pROC); library(mice); library(ggplot2)
})

DATA_FILE <- Sys.getenv("PANACEA_DATA", "5v_cleandf.RData")   # path to the downloaded source file
OUT_DIR   <- "revision_analyses/output"
CACHE_DIR <- "revision_analyses/cache"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE); dir.create(CACHE_DIR, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(OUT_DIR, f)

SEED_SPLIT <- 123; SEED_MICE <- 789
score_variables <- c("esi", "age", "meds_cardiovascular", "meds_analgesics",
                     "meds_psychotherapeuticdrugs", "n_admissions", "triage_vital_hr")
PUBLISHED_WEIGHTS <- c(esi = 3, age = 1, meds_cardiovascular = 5, meds_analgesics = 5,
                       meds_psychotherapeuticdrugs = 4, n_admissions = 2, triage_vital_hr = 2)
LOW_CUT <- 13; HIGH_CUT <- 21          # risk-band boundaries, derived on the derivation partition (script 01)

# ------------------------------------------------------------------------------
# 1. Dataset assembly and split (identical to Sections 1-2 of the main pipeline)
# ------------------------------------------------------------------------------
load(DATA_FILE)                                    # object `df`: 560,486 x 972
N_SOURCE <- nrow(df); N_DISPOSITION_MISSING <- sum(is.na(df$disposition))   # for the participant flow
df <- df[!is.na(df$disposition), ]
df <- df[, colSums(is.na(df)) != nrow(df)]
df$id <- seq_len(nrow(df))

binary_indicators <- df[, sapply(df, function(col) is.numeric(col) && max(col, na.rm = TRUE) == 1)] %>%
  mutate(across(everything(), as.factor))
demographics   <- df[, c("esi", "age", "gender", "maritalstatus", "employstatus", "disposition", "arrivalmode")]
triage_hx      <- df[, c("triage_vital_hr", "triage_vital_sbp", "triage_vital_dbp", "triage_vital_rr", "triage_vital_o2",
                         "triage_vital_o2_device", "triage_vital_temp", "previousdispo", "n_edvisits", "n_admissions", "n_surgeries")]
medications    <- df[, grep("^meds", names(df))]
chief_complaint <- df[, grep("^cc", names(df))]
past_history <- binary_indicators[, !(colnames(binary_indicators) %in% colnames(chief_complaint))]
past_history <- past_history[, !(grepl("last$", colnames(past_history)) | grepl("^meds", colnames(past_history)))]
past_history <- past_history[, !colnames(past_history) %in% c("bloodculture,routine_npos", "triage_vital_o2_device",
                                                              "o2_device_min", "o2_device_max", "o2_device_median")]
chief_complaint$disposition <- df$disposition; past_history$disposition <- df$disposition; medications$disposition <- df$disposition
analysis_data <- cbind(past_history, chief_complaint, medications, demographics, triage_hx)
analysis_data <- analysis_data[, !duplicated(as.list(analysis_data))]
analysis_data$id <- df$id
stopifnot(nrow(analysis_data) == 560486)

set.seed(SEED_SPLIT)
train_idx <- createDataPartition(analysis_data$disposition, p = 0.70, list = FALSE)
train_data_raw <- analysis_data[train_idx, ]
test_data_raw  <- analysis_data[-train_idx, ]
stopifnot(nrow(train_data_raw) == 392341, nrow(test_data_raw) == 168145)
y_train <- ifelse(train_data_raw$disposition == "Admit", 1, 0)
y_test  <- ifelse(test_data_raw$disposition  == "Admit", 1, 0)

# ------------------------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------------------------
impute_median_mode <- function(train, test, vars = NULL) {
  # Median/mode imputation with parameters estimated on `train` and applied to both.
  if (is.null(vars)) vars <- names(train)
  tr <- train; te <- test
  for (col in vars) {
    if (col %in% c("id", "disposition", "disposition_bin") || !(col %in% names(train))) next
    if (sum(is.na(train[[col]])) > 0 || sum(is.na(test[[col]])) > 0) {
      value <- if (is.numeric(train[[col]])) median(train[[col]], na.rm = TRUE) else names(sort(table(train[[col]]), decreasing = TRUE))[1]
      tr[[col]][is.na(train[[col]])] <- value; te[[col]][is.na(test[[col]])] <- value
    }
  }
  list(train = tr, test = te)
}

prepare_xgb_data <- function(data) {
  d <- data; d$disposition <- ifelse(d$disposition == "Admit", 1, 0); d <- d[, !names(d) %in% "id"]
  for (col in names(d)) if (is.factor(d[[col]]) || is.character(d[[col]])) d[[col]] <- as.numeric(as.factor(d[[col]]))
  d
}

calculate_score <- function(data, weights, cutoffs) {
  # PAnACEA-H: ESI (5 - ESI) x weight; age one weight unit per completed decade, capped at 9;
  # the three medication items scored when the count is > 0; prior admissions and heart rate at Youden cutoffs.
  s <- (5 - as.numeric(as.character(data$esi))) * weights[["esi"]]
  s <- s + pmin(floor(data$age / 10), 9) * weights[["age"]]
  for (med in c("meds_cardiovascular", "meds_analgesics", "meds_psychotherapeuticdrugs")) s <- s + ifelse(data[[med]] > 0, weights[[med]], 0)
  s + ifelse(data$n_admissions > cutoffs[["n_admissions"]], weights[["n_admissions"]], 0) +
      ifelse(data$triage_vital_hr > cutoffs[["triage_vital_hr"]], weights[["triage_vital_hr"]], 0)
}

find_youden_cutoff <- function(x, outcome) unname(coords(roc(outcome, x, quiet = TRUE), "best", ret = "threshold")$threshold)

fast_auc <- function(y, p) {                       # rank-based AUC, safe for very large n
  ok <- !is.na(p) & !is.na(y); y <- as.numeric(y[ok]); p <- p[ok]
  n1 <- as.double(sum(y == 1)); n0 <- as.double(length(y)) - n1
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

wilson_ci <- function(x, n, conf = 0.95) {
  if (n == 0) return(c(NA_real_, NA_real_))
  z <- qnorm(1 - (1 - conf) / 2); p <- x / n; d <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / d; hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(max(0, ctr - hw), min(1, ctr + hw))
}

calibration_metrics <- function(y, p, label = "") {
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10); lp <- qlogis(p)
  m_slope <- glm(y ~ lp, family = binomial); m_int <- glm(y ~ offset(lp), family = binomial)
  s_ci <- suppressMessages(confint.default(m_slope)[2, ]); i_ci <- suppressMessages(confint.default(m_int)[1, ])
  data.frame(Partition = label, N = length(y), Observed = sum(y), Expected = sum(p), OE_ratio = sum(y) / sum(p),
             Cal_intercept = unname(coef(m_int)[1]), Cal_int_low = i_ci[1], Cal_int_high = i_ci[2],
             Cal_slope = unname(coef(m_slope)[2]), Cal_slope_low = s_ci[1], Cal_slope_high = s_ci[2],
             Brier = mean((p - y)^2), Brier_scaled = 1 - mean((p - y)^2) / (mean(y) * (1 - mean(y))),
             AUC = as.numeric(auc(roc(y, p, quiet = TRUE))), stringsAsFactors = FALSE)
}

net_benefit <- function(y, p, t) { sel <- p >= t; sum(y[sel] == 1) / length(y) - (sum(y[sel] == 0) / length(y)) * (t / (1 - t)) }

# ------------------------------------------------------------------------------
# 3. Derivation of the score (multiple imputation, cutoffs, weights); cached
# ------------------------------------------------------------------------------
cache_file <- file.path(CACHE_DIR, "derivation.rds")
if (file.exists(cache_file)) {
  D <- readRDS(cache_file)
} else {
  cat("Deriving the score on the derivation partition (multiple imputation; this takes a while)...\n")
  tr_scores <- train_data_raw[, c(score_variables, "disposition", "id")]
  fit <- mice(tr_scores[, score_variables], m = 5, maxit = 10, method = "pmm", seed = SEED_MICE, printFlag = FALSE)
  train_imputed <- complete(fit, action = 1)                    # first completed dataset (not pooled), as in the main pipeline
  train_imputed$disposition_bin <- y_train; train_imputed$id <- tr_scores$id
  cutoffs <- list(n_admissions = find_youden_cutoff(train_imputed$n_admissions, y_train),
                  triage_vital_hr = find_youden_cutoff(train_imputed$triage_vital_hr, y_train))
  weights <- setNames(as.list(rep(1, length(score_variables))), score_variables)
  for (v in score_variables) {
    best_w <- 1; best_a <- -Inf
    for (w in 1:5) { cand <- weights; cand[[v]] <- w
      a <- fast_auc(y_train, calculate_score(train_imputed, cand, cutoffs)); if (a > best_a) { best_a <- a; best_w <- w } }
    weights[[v]] <- best_w
  }
  imputation_params <- lapply(score_variables, function(v) if (is.numeric(train_imputed[[v]])) median(train_imputed[[v]], na.rm = TRUE) else names(sort(table(train_imputed[[v]]), decreasing = TRUE))[1])
  names(imputation_params) <- score_variables
  D <- list(train_imputed = train_imputed, cutoffs = cutoffs, weights = weights, imputation_params = imputation_params)
  saveRDS(D, cache_file)
}
train_imputed <- D$train_imputed; cutoffs <- D$cutoffs; weights <- D$weights; imputation_params <- D$imputation_params
stopifnot(all(unlist(weights[score_variables]) == PUBLISHED_WEIGHTS[score_variables]))    # the derivation reproduces Table 2
train_imputed$score <- calculate_score(train_imputed, weights, cutoffs)

# held-out partition: median/mode values of the derivation partition (no multiple imputation)
test_final <- test_data_raw[, c(score_variables, "disposition", "id")]
for (v in score_variables) if (anyNA(test_final[[v]])) test_final[[v]][is.na(test_final[[v]])] <- imputation_params[[v]]
test_final$disposition_bin <- y_test
test_final$score <- calculate_score(test_final, weights, cutoffs)
stopifnot(abs(fast_auc(y_test, test_final$score) - 0.8797) < 5e-4)                        # Table 3 / Abstract
cat(sprintf("Setup complete: derivation n = %d, held-out n = %d, held-out AUC = %.4f\n",
            nrow(train_imputed), nrow(test_final), fast_auc(y_test, test_final$score)))
