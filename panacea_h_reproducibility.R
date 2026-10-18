# ==============================================================================
# PAnACEA-H: development and internal validation of a seven-item score for
# predicting hospital admission at Emergency Department triage
#
# This script reproduces, from the publicly available dataset, every model-
# development and validation result reported in the manuscript:
#   - Feature ranking by information gain (Figure 2; Supplementary Table 4)
#   - Pooled 5-fold cross-validated performance of the XGBoost ranking model
#     (Supplementary Table 5)
#   - The PAnACEA-H score itself: item weights and cutoffs (Table 2)
#   - Internal validation by 10-fold cross-validation on the training set
#   - Final validation on the held-out test set, including the ROC curve
#     (Figure 3), the threshold-level performance table (Table 3), and
#     bootstrap confidence intervals (Supplementary Table 6, Supplementary
#     Figure 1)
#   - Decision curve analysis (Supplementary Figure 2)
#
# NOT covered by this script: Table 1 (descriptive baseline characteristics
# by disposition). That table is purely descriptive, does not feed into the
# score's derivation or validation, and is not reproduced here.
#
# --------------------------------------------------------------------------
# DATA
# --------------------------------------------------------------------------
# Source study:
#   Hong WS, Haimovich AD, Taylor RA (2018). Predicting hospital admission
#   at emergency department triage using machine learning. PLoS ONE 13(7):
#   e0201016. https://doi.org/10.1371/journal.pone.0201016
#
# Dataset (de-identified, 560,486 visits x 972 variables), publicly
# available at:
#   https://github.com/yaleemmlc/admissionprediction
#   -> Results/5v_cleandf.RData
#
# Download that file and point DATA_FILE (below) to its local path before
# running this script.
#
# --------------------------------------------------------------------------
# SOFTWARE
# --------------------------------------------------------------------------
# R (this script was run under R 4.2.3). Required packages:
#   dplyr, tidyr, caret, xgboost, pROC, mice, ggplot2, dcurves
#
# Note on the decision-curve package: the manuscript's own analysis used the
# `rmda` package, which has since been withdrawn from CRAN. `dcurves` is the
# actively maintained replacement and is used here instead; unlike `rmda`,
# `dcurves` does not fit its own logistic model internally, so this script
# fits one explicitly (see Section 9) to obtain predicted probabilities
# before calling dcurves::dca().
#
# --------------------------------------------------------------------------
# REPRODUCIBILITY NOTE ON XGBOOST HYPERPARAMETERS (max_depth)
# --------------------------------------------------------------------------
# The ranking model's hyperparameters (max_depth = 3, eta = 0.1, 100 rounds)
# are fixed below to the values stated in the Methods, not selected by the
# grid search in Section 4b. That grid search is included for transparency:
# run at this sample size, it identifies max_depth = 8 as achieving higher
# cross-validated AUC (~0.92) than max_depth = 3 (~0.914, the value that
# matches the manuscript's reported 91.4%). The published pipeline used the
# fixed value, not the grid search's own recommendation; this discrepancy is
# disclosed in the response to reviewers accompanying this submission.
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(caret)
library(xgboost)
library(pROC)
library(mice)
library(ggplot2)
library(dcurves)

DATA_FILE <- "5v_cleandf.RData"   # <- edit to the local path of the downloaded file
OUTPUT_DIR <- "reproducibility_output"
dir.create(OUTPUT_DIR, showWarnings = FALSE)
out <- function(filename) file.path(OUTPUT_DIR, filename)

SEED_SUBSET_MERGE   <- NULL   # no subsampling is used anywhere in this script
SEED_TRAIN_TEST_SPLIT <- 123
SEED_XGB_FOLDS        <- 456
SEED_MICE             <- 789
SEED_SCORE_CV_FOLDS   <- 111
SEED_BOOTSTRAP        <- 2024

# ------------------------------------------------------------------------------
# 1. Load and assemble the analysis dataset (547 variables, 560,486 visits)
# ------------------------------------------------------------------------------
cat("=== 1. LOADING AND ASSEMBLING THE DATASET ===\n")

load(DATA_FILE)   # loads an object named `df`: 560,486 x 972

df <- df[!is.na(df$disposition), ]
df <- df[, colSums(is.na(df)) != nrow(df)]
df$id <- seq_len(nrow(df))

binary_indicators <- df[, sapply(df, function(col) is.numeric(col) && max(col, na.rm = TRUE) == 1)] %>%
  mutate(across(everything(), as.factor))

demographics_cols <- c("esi", "age", "gender", "maritalstatus", "employstatus",
                       "disposition", "arrivalmode")
triage_hx_cols <- c("triage_vital_hr", "triage_vital_sbp", "triage_vital_dbp",
                    "triage_vital_rr", "triage_vital_o2", "triage_vital_o2_device",
                    "triage_vital_temp", "previousdispo", "n_edvisits",
                    "n_admissions", "n_surgeries")

demographics   <- df[, demographics_cols]
triage_hx      <- df[, triage_hx_cols]
medications    <- df[, grep("^meds", names(df))]
chief_complaint <- df[, grep("^cc", names(df))]

past_history <- binary_indicators[, !(colnames(binary_indicators) %in% colnames(chief_complaint))]
past_history <- past_history[, !(grepl("last$", colnames(past_history)) | grepl("^meds", colnames(past_history)))]
past_history <- past_history[, !colnames(past_history) %in% c(
  "bloodculture,routine_npos", "triage_vital_o2_device",
  "o2_device_min", "o2_device_max", "o2_device_median")]

chief_complaint$disposition <- df$disposition
past_history$disposition    <- df$disposition
medications$disposition     <- df$disposition

analysis_data <- cbind(past_history, chief_complaint, medications, demographics, triage_hx)
analysis_data <- analysis_data[, !duplicated(as.list(analysis_data))]
analysis_data$id <- df$id

cat(sprintf("Assembled dataset: %d visits x %d variables\n",
            nrow(analysis_data), ncol(analysis_data)))
stopifnot(nrow(analysis_data) == 560486)

cat(sprintf("Admissions: %d (%.1f%%)\n",
            sum(analysis_data$disposition == "Admit"),
            100 * mean(analysis_data$disposition == "Admit")))

# ------------------------------------------------------------------------------
# 2. Derivation / validation split (70% / 30%)
# ------------------------------------------------------------------------------
cat("\n=== 2. DERIVATION / VALIDATION SPLIT ===\n")

set.seed(SEED_TRAIN_TEST_SPLIT)
train_idx <- createDataPartition(analysis_data$disposition, p = 0.70, list = FALSE)
train_data_raw <- analysis_data[train_idx, ]
test_data_raw  <- analysis_data[-train_idx, ]

cat(sprintf("Derivation (training) set: %d visits (%.1f%% admitted)\n",
            nrow(train_data_raw), 100 * mean(train_data_raw$disposition == "Admit")))
cat(sprintf("Validation (test) set:     %d visits (%.1f%% admitted)\n",
            nrow(test_data_raw), 100 * mean(test_data_raw$disposition == "Admit")))
stopifnot(nrow(train_data_raw) == 392341, nrow(test_data_raw) == 168145)

# ------------------------------------------------------------------------------
# 3. Helper functions
# ------------------------------------------------------------------------------
impute_median_mode <- function(train, test, vars = NULL) {
  # Median/mode imputation, parameters estimated on `train` and applied to
  # both `train` and `test` (no information from `test` is used).
  if (is.null(vars)) vars <- names(train)
  train_imputed <- train
  test_imputed  <- test
  params <- list()
  for (col in vars) {
    if (col %in% c("id", "disposition", "disposition_bin")) next
    if (!(col %in% names(train))) next
    if (sum(is.na(train[[col]])) > 0 || sum(is.na(test[[col]])) > 0) {
      if (is.numeric(train[[col]])) {
        value <- median(train[[col]], na.rm = TRUE)
      } else {
        value <- names(sort(table(train[[col]]), decreasing = TRUE))[1]
      }
      params[[col]] <- value
      train_imputed[[col]][is.na(train[[col]])] <- value
      test_imputed[[col]][is.na(test[[col]])]   <- value
    }
  }
  list(train = train_imputed, test = test_imputed, params = params)
}

prepare_xgb_data <- function(data) {
  # Numeric label + numeric-coded predictor matrix for xgboost.
  data_xgb <- data
  data_xgb$disposition <- ifelse(data_xgb$disposition == "Admit", 1, 0)
  data_xgb <- data_xgb[, !names(data_xgb) %in% "id"]
  for (col in names(data_xgb)) {
    if (is.factor(data_xgb[[col]]) || is.character(data_xgb[[col]])) {
      data_xgb[[col]] <- as.numeric(as.factor(data_xgb[[col]]))
    }
  }
  data_xgb
}

calculate_score <- function(data, weights, cutoffs) {
  # PAnACEA-H integer score.
  # ESI: multiplicative, base weight x ordinal distance from ESI 5.
  # Age: cumulative, one weight unit per completed decade, capped at 9
  #      (i.e. patients aged 90+ receive the maximum age contribution).
  # Cardiovascular / analgesic / psychiatric medications: dichotomous, full
  #      weight awarded when the corresponding count is > 0.
  # Prior-year admissions and triage heart rate: dichotomous at a
  #      Youden-optimal cutoff (see find_youden_cutoff()).
  score <- (5 - as.numeric(as.character(data$esi))) * weights[["esi"]]
  age_decade <- pmin(floor(data$age / 10), 9)
  score <- score + age_decade * weights[["age"]]
  for (med in c("meds_cardiovascular", "meds_analgesics", "meds_psychotherapeuticdrugs")) {
    score <- score + ifelse(data[[med]] > 0, weights[[med]], 0)
  }
  score <- score + ifelse(data$n_admissions > cutoffs[["n_admissions"]],
                          weights[["n_admissions"]], 0)
  score <- score + ifelse(data$triage_vital_hr > cutoffs[["triage_vital_hr"]],
                          weights[["triage_vital_hr"]], 0)
  score
}

find_youden_cutoff <- function(x, outcome) {
  r <- roc(outcome, x, quiet = TRUE)
  unname(coords(r, "best", ret = "threshold")$threshold)
}

evaluate_score_auc <- function(weights, data, cutoffs) {
  s <- calculate_score(data, weights, cutoffs)
  if (length(unique(s)) > 1) as.numeric(auc(roc(data$disposition_bin, s, quiet = TRUE))) else 0.5
}

calculate_all_metrics <- function(actual, predicted, threshold = NULL) {
  # `predicted` may be a probability or, for the score, the integer value
  # itself: coords()/roc() only use its ranking, so both work identically.
  r <- roc(actual, predicted, quiet = TRUE)
  auc_val <- as.numeric(auc(r))
  if (is.null(threshold)) {
    threshold <- coords(r, "best", ret = "threshold")$threshold
  }
  pred_class <- ifelse(predicted >= threshold, 1, 0)
  tp <- sum(pred_class == 1 & actual == 1); tn <- sum(pred_class == 0 & actual == 0)
  fp <- sum(pred_class == 1 & actual == 0); fn <- sum(pred_class == 0 & actual == 1)
  list(AUC = auc_val, Threshold = threshold,
       Accuracy = (tp + tn) / (tp + tn + fp + fn),
       Sensitivity = tp / (tp + fn), Specificity = tn / (tn + fp),
       PPV = tp / (tp + fp), NPV = tn / (tn + fn))
}

# ==============================================================================
# 4. Feature ranking: XGBoost, 5-fold cross-validation
# ==============================================================================
cat("\n=== 4. FEATURE RANKING (XGBoost, 5-fold CV) ===\n")

k_xgb <- 5
set.seed(SEED_XGB_FOLDS)
xgb_folds <- createFolds(train_data_raw$disposition, k = k_xgb, list = TRUE, returnTrain = FALSE)

# --- 4a. Fixed hyperparameters used for the published ranking -----------------
XGB_MAX_DEPTH <- 3
XGB_ETA <- 0.1
XGB_NROUNDS <- 100

importance_by_fold <- list()
auc_by_fold <- numeric(k_xgb)

for (i in seq_len(k_xgb)) {
  test_idx <- xgb_folds[[i]]
  fold_train_raw <- train_data_raw[-test_idx, ]
  fold_test_raw  <- train_data_raw[test_idx, ]

  imputed <- impute_median_mode(fold_train_raw, fold_test_raw)
  fold_train <- prepare_xgb_data(imputed$train)
  fold_test  <- prepare_xgb_data(imputed$test)

  x_train <- as.matrix(fold_train[, setdiff(names(fold_train), "disposition")])
  x_test  <- as.matrix(fold_test[,  setdiff(names(fold_test),  "disposition")])
  dtrain <- xgb.DMatrix(x_train, label = fold_train$disposition)
  dtest  <- xgb.DMatrix(x_test,  label = fold_test$disposition)

  fold_model <- xgb.train(
    params = list(objective = "binary:logistic", eval_metric = "auc",
                  max_depth = XGB_MAX_DEPTH, eta = XGB_ETA),
    data = dtrain, nrounds = XGB_NROUNDS, verbose = 0)

  pred <- predict(fold_model, dtest)
  auc_by_fold[i] <- as.numeric(auc(roc(fold_test$disposition, pred, quiet = TRUE)))

  imp <- xgb.importance(feature_names = colnames(x_train), model = fold_model)
  imp$Fold <- i
  importance_by_fold[[i]] <- imp

  cat(sprintf("  Fold %d/%d: AUC = %.4f\n", i, k_xgb, auc_by_fold[i]))
}

cat(sprintf("\nPooled 5-fold cross-validated AUC: %.4f (SD %.4f)\n",
            mean(auc_by_fold), sd(auc_by_fold)))
cat("Manuscript-reported value: 0.9142 (Results: \"a robust pooled AUC of 91.4%\")\n")

feature_ranking <- bind_rows(importance_by_fold) %>%
  group_by(Feature) %>%
  summarise(Frequency = n(), Mean_Gain = mean(Gain, na.rm = TRUE),
            SD_Gain = sd(Gain, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(Mean_Gain))

cat("\nTop 15 variables by mean information gain (Figure 2 / Supplementary Table 4):\n")
print(head(feature_ranking, 15), n = 15)

write.csv(data.frame(Fold = c(paste0("Fold_", seq_len(k_xgb)), "Pooled_Mean"),
                     AUC = c(auc_by_fold, mean(auc_by_fold))),
         out("supplementary_table5_xgboost_cv_auc.csv"), row.names = FALSE)
write.csv(feature_ranking, out("supplementary_table4_feature_ranking.csv"), row.names = FALSE)

p_importance <- ggplot(head(feature_ranking, 15),
                       aes(x = reorder(Feature, Mean_Gain), y = Mean_Gain)) +
  geom_point(size = 3, colour = "steelblue") +
  geom_errorbar(aes(ymin = Mean_Gain - SD_Gain, ymax = Mean_Gain + SD_Gain), width = 0.2) +
  coord_flip() +
  labs(title = "Feature importance (5-fold cross-validation)",
       subtitle = sprintf("XGBoost pooled AUC: %.3f (SD %.3f)", mean(auc_by_fold), sd(auc_by_fold)),
       x = NULL, y = "Mean information gain (+/- SD)") +
  theme_minimal()
ggsave(out("figure2_feature_importance.jpeg"), p_importance, dpi = 300, bg = "white",
       width = 9, height = 7)

# --- 4b. Grid search over max_depth (diagnostic only; see header note) --------
cat("\n--- 4b. Diagnostic: max_depth grid search (not used for the published model) ---\n")
full_imputed <- impute_median_mode(train_data_raw, train_data_raw)
full_xgb <- prepare_xgb_data(full_imputed$train)
x_full <- as.matrix(full_xgb[, setdiff(names(full_xgb), "disposition")])
dfull  <- xgb.DMatrix(x_full, label = full_xgb$disposition)

depth_grid <- c(2, 3, 4, 5, 6, 8)
depth_grid_auc <- sapply(depth_grid, function(d) {
  cv <- xgb.cv(data = dfull, nrounds = XGB_NROUNDS, folds = xgb_folds,
              params = list(objective = "binary:logistic", eval_metric = "auc",
                            max_depth = d, eta = XGB_ETA), verbose = 0)
  cv$evaluation_log$test_auc_mean[XGB_NROUNDS]
})
depth_grid_result <- data.frame(max_depth = depth_grid, cv_auc = depth_grid_auc)
print(depth_grid_result, row.names = FALSE)
cat(sprintf("Grid maximum: max_depth = %d (CV AUC %.4f). Published model uses max_depth = %d.\n",
            depth_grid_result$max_depth[which.max(depth_grid_result$cv_auc)],
            max(depth_grid_result$cv_auc), XGB_MAX_DEPTH))
write.csv(depth_grid_result, out("diagnostic_max_depth_grid_search.csv"), row.names = FALSE)

# ==============================================================================
# 5. Score derivation
# ==============================================================================
cat("\n\n=== 5. SCORE DERIVATION ===\n")

score_variables <- c("esi", "age", "meds_cardiovascular", "meds_analgesics",
                     "meds_psychotherapeuticdrugs", "n_admissions", "triage_vital_hr")
cat("Selected predictors (clinical relevance and cross-system reproducibility,\n",
    "applied to the top-15 ranking above):", paste(score_variables, collapse = ", "), "\n")

# --- 5a. Multiple imputation of the seven score variables on the full
#         training set (MICE, m = 5, predictive mean matching) ---------------
train_score_raw <- train_data_raw[, c(score_variables, "disposition", "id")]

cat("\nRunning multiple imputation (mice, m = 5, maxit = 10, pmm)...\n")
mice_fit <- mice(train_score_raw[, score_variables], m = 5, maxit = 10,
                 method = "pmm", seed = SEED_MICE, printFlag = FALSE)
train_imputed <- complete(mice_fit, action = 1)
train_imputed$disposition     <- train_score_raw$disposition
train_imputed$id              <- train_score_raw$id
train_imputed$disposition_bin <- ifelse(train_imputed$disposition == "Admit", 1, 0)

# --- 5b. Cutoffs (Youden's index on the imputed training set, fixed before
#         weight optimisation) -------------------------------------------------
cutoffs <- list(
  n_admissions    = find_youden_cutoff(train_imputed$n_admissions, train_imputed$disposition_bin),
  triage_vital_hr = find_youden_cutoff(train_imputed$triage_vital_hr, train_imputed$disposition_bin)
)
cat(sprintf("\nCutoffs (Youden's index, training set): n_admissions > %.2f, triage_vital_hr > %.2f\n",
            cutoffs$n_admissions, cutoffs$triage_vital_hr))

# --- 5c. Greedy coordinate-ascent weight optimisation --------------------------
# Variables are visited once, in the order given by `score_variables`; for
# each, integer weights 1-5 are tried holding all other weights fixed, the
# AUC-maximising weight is kept, and the algorithm moves on (no second pass).
weights <- setNames(as.list(rep(1, length(score_variables))), score_variables)
for (v in score_variables) {
  best_w <- 1; best_auc <- -Inf
  for (w in 1:5) {
    candidate <- weights; candidate[[v]] <- w
    a <- evaluate_score_auc(candidate, train_imputed, cutoffs)
    if (a > best_auc) { best_auc <- a; best_w <- w }
  }
  weights[[v]] <- best_w
}

cat("\nDerived weights (Table 2):\n")
print(unlist(weights))
cat(sprintf("Apparent (training) AUC at these weights: %.4f\n",
            evaluate_score_auc(weights, train_imputed, cutoffs)))
cat("Manuscript-published weights: ESI 3, age 1/decade, cardiovascular 5,",
    "analgesics 5, psychiatric 4, prior admissions 2, heart rate 2\n")
cat(sprintf("Maximum attainable score: %d points\n",
            12 + 9 + weights$meds_cardiovascular + weights$meds_analgesics +
              weights$meds_psychotherapeuticdrugs + weights$n_admissions +
              weights$triage_vital_hr))

write.csv(data.frame(Item = names(weights), Weight = unlist(weights)),
         out("table2_score_weights.csv"), row.names = FALSE)

# ==============================================================================
# 6. Internal validation: 10-fold cross-validation on the training set
# ==============================================================================
cat("\n\n=== 6. INTERNAL VALIDATION (10-fold CV, weights and cutoffs held fixed) ===\n")
# MICE is repeated inside each fold (imputation uncertainty); the item
# weights and cutoffs derived in Section 5 are NOT re-estimated per fold.
# As discussed in the manuscript and the accompanying response to reviewers,
# this is a robustness check with respect to imputation and sampling
# variability, not an independent validation of the derivation procedure.

k_score <- 10
set.seed(SEED_SCORE_CV_FOLDS)
score_folds <- createFolds(train_data_raw$disposition, k = k_score, list = TRUE, returnTrain = FALSE)

cv_metrics <- vector("list", k_score)
for (i in seq_len(k_score)) {
  test_idx <- score_folds[[i]]
  fold_train_raw <- train_data_raw[-test_idx, c(score_variables, "disposition", "id")]
  fold_test_raw  <- train_data_raw[test_idx,  c(score_variables, "disposition", "id")]

  fold_mice <- mice(fold_train_raw[, score_variables], m = 5, maxit = 10,
                    method = "pmm", seed = SEED_MICE + i, printFlag = FALSE)
  fold_train_imp <- complete(fold_mice, action = 1)
  fold_train_imp$disposition_bin <- ifelse(fold_train_raw$disposition == "Admit", 1, 0)

  fold_test_imp <- fold_test_raw
  for (v in score_variables) {
    if (anyNA(fold_test_raw[[v]])) {
      value <- if (is.numeric(fold_train_imp[[v]])) median(fold_train_imp[[v]], na.rm = TRUE)
      else names(sort(table(fold_train_imp[[v]]), decreasing = TRUE))[1]
      fold_test_imp[[v]][is.na(fold_test_raw[[v]])] <- value
    }
  }
  fold_test_imp$disposition_bin <- ifelse(fold_test_raw$disposition == "Admit", 1, 0)

  fold_score <- calculate_score(fold_test_imp, weights, cutoffs)
  cv_metrics[[i]] <- as.data.frame(c(calculate_all_metrics(fold_test_imp$disposition_bin, fold_score),
                                     Fold = i))
  cat(sprintf("  Fold %d/%d done.\n", i, k_score))
}
cv_results <- bind_rows(cv_metrics)

cv_summary <- cv_results %>%
  summarise(across(c(AUC, Accuracy, Sensitivity, Specificity, PPV, NPV),
                   list(mean = mean, sd = sd)))
cat("\n10-fold CV performance (mean +/- SD):\n")
print(cv_summary)
cat("Manuscript-reported: AUC 88.1% (95% CI 88.0-88.2%), accuracy 80.8%,",
    "sensitivity 74.7%, specificity 83.4%, PPV 65.6%, NPV 88.6%\n")

write.csv(cv_results, out("internal_cv_by_fold.csv"), row.names = FALSE)
write.csv(cv_summary, out("internal_cv_summary.csv"), row.names = FALSE)

# ==============================================================================
# 7. Final validation on the held-out test set
# ==============================================================================
cat("\n\n=== 7. HELD-OUT TEST SET VALIDATION ===\n")
# Missing values are imputed with the median/mode of the (imputed) training
# set -- no information from the test set is used, and MICE is not applied
# here (as in the manuscript, for a bedside-realistic imputation strategy).

imputation_params <- lapply(score_variables, function(v) {
  if (is.numeric(train_imputed[[v]])) median(train_imputed[[v]], na.rm = TRUE)
  else names(sort(table(train_imputed[[v]]), decreasing = TRUE))[1]
})
names(imputation_params) <- score_variables

test_final <- test_data_raw[, c(score_variables, "disposition", "id")]
for (v in score_variables) {
  if (anyNA(test_final[[v]])) test_final[[v]][is.na(test_final[[v]])] <- imputation_params[[v]]
}
test_final$disposition_bin <- ifelse(test_final$disposition == "Admit", 1, 0)
test_final$score <- calculate_score(test_final, weights, cutoffs)

test_metrics <- calculate_all_metrics(test_final$disposition_bin, test_final$score)
roc_test <- roc(test_final$disposition_bin, test_final$score, quiet = TRUE)
ci_test  <- ci.auc(roc_test)

cat(sprintf("\nTest-set AUC: %.4f (95%% CI %.4f-%.4f)\n", test_metrics$AUC, ci_test[1], ci_test[3]))
cat(sprintf("At the best cutoff (score >= %.0f): accuracy %.4f, sensitivity %.4f, specificity %.4f, PPV %.4f, NPV %.4f\n",
            test_metrics$Threshold, test_metrics$Accuracy, test_metrics$Sensitivity,
            test_metrics$Specificity, test_metrics$PPV, test_metrics$NPV))
cat("Manuscript-reported: AUC 88.0% (95% CI 87.8-88.1%)\n")

p_roc <- ggroc(roc_test, colour = "darkblue", linewidth = 1.1) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", colour = "grey50") +
  labs(title = "ROC curve, held-out test set",
       subtitle = sprintf("AUC = %.3f (95%% CI %.3f-%.3f)", test_metrics$AUC, ci_test[1], ci_test[3]),
       x = "Specificity", y = "Sensitivity") +
  theme_minimal()
ggsave(out("figure3_roc_test_set.jpeg"), p_roc, dpi = 300, bg = "white", width = 7, height = 6)

# --- Table 3: performance at every attainable score threshold -----------------
thresholds <- sort(unique(test_final$score))
threshold_table <- lapply(thresholds, function(t) {
  pred <- as.integer(test_final$score >= t)
  y <- test_final$disposition_bin
  tp <- sum(pred == 1 & y == 1); tn <- sum(pred == 0 & y == 0)
  fp <- sum(pred == 1 & y == 0); fn <- sum(pred == 0 & y == 1)
  data.frame(Threshold = t, Sensitivity = tp / (tp + fn), Specificity = tn / (tn + fp),
            PPV = tp / (tp + fp), NPV = tn / (tn + fn),
            Accuracy = (tp + tn) / nrow(test_final))
})
threshold_table <- bind_rows(threshold_table)
cat("\nTable 3 (threshold-level performance, held-out test set):\n")
print(threshold_table, digits = 3)
write.csv(threshold_table, out("table3_threshold_performance.csv"), row.names = FALSE)

low_band  <- threshold_table$Threshold[max(which(threshold_table$NPV >= 0.90))]
high_band <- threshold_table$Threshold[min(which(threshold_table$PPV >= 0.90))]
cat(sprintf("\nRisk bands: Low 0-%d | Intermediate %d-%d | High >= %d\n",
            low_band, low_band + 1, high_band - 1, high_band))
cat("Manuscript-published bands: Low 0-14 | Intermediate 15-20 | High >= 21\n")

# ==============================================================================
# 8. Bootstrap confidence intervals on the held-out test set
# ==============================================================================
cat("\n\n=== 8. BOOTSTRAP VALIDATION (1000 replicates) ===\n")

set.seed(SEED_BOOTSTRAP)
n_boot <- 1000
n_test <- nrow(test_final)
boot_metrics <- matrix(NA_real_, n_boot, 6,
                       dimnames = list(NULL, c("AUC", "Accuracy", "Sensitivity",
                                               "Specificity", "PPV", "NPV")))
for (b in seq_len(n_boot)) {
  idx <- sample.int(n_test, n_test, replace = TRUE)
  m <- calculate_all_metrics(test_final$disposition_bin[idx], test_final$score[idx],
                             threshold = test_metrics$Threshold)
  boot_metrics[b, ] <- unlist(m[c("AUC", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")])
}
boot_summary <- data.frame(
  Metric = colnames(boot_metrics),
  Estimate = unlist(test_metrics[c("AUC", "Accuracy", "Sensitivity", "Specificity", "PPV", "NPV")]),
  CI_low  = apply(boot_metrics, 2, quantile, 0.025),
  CI_high = apply(boot_metrics, 2, quantile, 0.975))
cat("Bootstrap 95% CIs (Supplementary Table 6):\n")
print(boot_summary, row.names = FALSE, digits = 4)
write.csv(boot_summary, out("supplementary_table6_bootstrap_ci.csv"), row.names = FALSE)
write.csv(as.data.frame(boot_metrics), out("supplementary_figure1_bootstrap_distributions.csv"),
         row.names = FALSE)

# ==============================================================================
# 9. Decision curve analysis
# ==============================================================================
cat("\n\n=== 9. DECISION CURVE ANALYSIS ===\n")
# dcurves requires predicted probabilities rather than the raw integer
# score; a single-predictor logistic model converts one into the other
# (see header note on the rmda -> dcurves substitution).
score_to_prob <- glm(disposition_bin ~ score, data = test_final, family = binomial)
test_final$predicted_prob <- predict(score_to_prob, type = "response")

dca_fit <- dcurves::dca(disposition_bin ~ predicted_prob, data = test_final,
                        thresholds = seq(0.01, 0.99, by = 0.01),
                        label = list(predicted_prob = "PAnACEA-H score"))
p_dca <- plot(dca_fit, smooth = FALSE) +
  labs(title = "Decision curve analysis, held-out test set", x = "Threshold probability") +
  theme_minimal()
ggsave(out("supplementary_figure2_decision_curve.jpeg"), p_dca, dpi = 300, bg = "white",
      width = 8, height = 6)
write.csv(as.data.frame(dca_fit$dca), out("supplementary_figure2_decision_curve_data.csv"),
         row.names = FALSE)

# ==============================================================================
# 10. Summary
# ==============================================================================
cat("\n\n========== SUMMARY ==========\n")
cat(sprintf("Dataset: %d visits (%d derivation, %d validation)\n",
            nrow(analysis_data), nrow(train_data_raw), nrow(test_data_raw)))
cat(sprintf("XGBoost ranking model, pooled 5-fold CV AUC: %.4f\n", mean(auc_by_fold)))
cat("Score weights:\n"); print(unlist(weights))
cat(sprintf("Score, 10-fold CV AUC: %.4f\n", cv_summary$AUC_mean))
cat(sprintf("Score, held-out test AUC: %.4f (95%% CI %.4f-%.4f)\n",
            test_metrics$AUC, ci_test[1], ci_test[3]))
cat(sprintf("\nAll output files written to ./%s/\n", OUTPUT_DIR))

save(weights, cutoffs, imputation_params, feature_ranking, auc_by_fold,
    cv_summary, test_metrics, ci_test, threshold_table, boot_summary,
    file = out("panacea_h_results.RData"))
