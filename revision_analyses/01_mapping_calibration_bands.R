# ==============================================================================
# 01. Score-to-probability mapping, calibration, risk bands, decision curve
#
# Reproduces (manuscript numbering):
#   Supplementary Table 17  score -> predicted probability (derivation partition only)
#   Supplementary Table 18  threshold-level metrics, derivation partition
#   Table 3                 threshold-level metrics, validation partition
#   Supplementary Figure 3  calibration (data behind the plot)
#   Supplementary Figure 4  decision curve (score, ESI alone, Admit all / none)
#   risk-band boundaries and the performance of the bands (Results; Supplementary Table 25)
#
# Principle: everything that defines the tool (the mapping, the band boundaries,
# the operating point) is estimated on the derivation partition only and then
# frozen before the held-out partition is touched.
#
# Run from the repository root:  Rscript revision_analyses/01_mapping_calibration_bands.R
# ==============================================================================
source("revision_analyses/R/setup.R")

# ---- 1. Empirical score-to-probability mapping (derivation partition only) ----
POOL_LOW <- 4; POOL_HIGH <- 33          # scores <= 4 and >= 33 share one estimate (sparse cells; keeps the map monotone)
emp_map <- train_imputed %>%
  mutate(score_grp = pmin(pmax(score, POOL_LOW), POOL_HIGH)) %>%
  group_by(score_grp) %>%
  summarise(n_train = n(), n_admit = sum(disposition_bin), p_map = mean(disposition_bin), .groups = "drop") %>%
  rowwise() %>% mutate(ci_low = wilson_ci(n_admit, n_train)[1], ci_high = wilson_ci(n_admit, n_train)[2]) %>% ungroup()
stopifnot(all(diff(emp_map$p_map) >= 0))                                  # monotone over the whole range
apply_map <- function(score) emp_map$p_map[match(pmin(pmax(score, POOL_LOW), POOL_HIGH), emp_map$score_grp)]

# logistic approximation, reported next to the empirical map
fit_logit <- glm(disposition_bin ~ score, data = train_imputed, family = binomial)
alpha <- unname(coef(fit_logit)[1]); beta <- unname(coef(fit_logit)[2])

# ---- 2. Threshold-level metrics, both partitions; risk bands from the derivation partition ----
threshold_table <- function(data, label) {
  bind_rows(lapply(sort(unique(data$score)), function(t) {
    pred <- as.integer(data$score >= t); y <- data$disposition_bin            # a threshold means: score greater than or equal to t
    tp <- sum(pred == 1 & y == 1); tn <- sum(pred == 0 & y == 0); fp <- sum(pred == 1 & y == 0); fn <- sum(pred == 0 & y == 1)
    data.frame(Threshold = t, Sensitivity = tp / (tp + fn), Specificity = tn / (tn + fp),
               PPV = ifelse(tp + fp > 0, tp / (tp + fp), NA), NPV = ifelse(tn + fn > 0, tn / (tn + fn), NA),
               Accuracy = (tp + tn) / nrow(data), Youden = tp / (tp + fn) + tn / (tn + fp) - 1)
  })) %>% mutate(Partition = label)
}
thr_train <- threshold_table(train_imputed, "Training"); thr_test <- threshold_table(test_final, "Test")

# pre-specified rules, applied to the derivation partition only:
#   low boundary  = highest score s such that the NPV for patients scoring <= s is at least 0.90
#   high boundary = lowest score t such that the PPV for patients scoring >= t is at least 0.90
npv_by_upper <- thr_train %>% mutate(upper = Threshold - 1) %>% filter(!is.na(NPV), upper >= min(train_imputed$score))
low_cut  <- max(npv_by_upper$upper[npv_by_upper$NPV >= 0.90])
high_cut <- min(thr_train$Threshold[!is.na(thr_train$PPV) & thr_train$PPV >= 0.90])
op_point <- thr_train$Threshold[which.max(thr_train$Youden)]               # Youden-optimal cutoff (reported for comparability only)
cat(sprintf("Risk bands: Low 0-%d | Intermediate %d-%d | High >= %d | Youden-optimal cutoff: score >= %d\n", low_cut, low_cut + 1, high_cut - 1, high_cut, op_point))
stopifnot(low_cut == LOW_CUT, high_cut == HIGH_CUT, op_point == 15)

publish_map <- data.frame(score = 0:39) %>%
  mutate(score_grp = pmin(pmax(score, POOL_LOW), POOL_HIGH)) %>% left_join(emp_map, by = "score_grp") %>%
  mutate(p_logistic_approx = plogis(alpha + beta * score),
         risk_band = case_when(score <= low_cut ~ "Low", score < high_cut ~ "Intermediate", TRUE ~ "High")) %>%
  select(score, risk_band, n_train, n_admit, p_predicted = p_map, ci_low, ci_high, p_logistic_approx)
write.csv(publish_map, out("supp_table17_score_to_probability.csv"), row.names = FALSE)
write.csv(bind_rows(thr_train, thr_test), out("supp_table18_table3_threshold_metrics_both_partitions.csv"), row.names = FALSE)

# ---- 3. Calibration -------------------------------------------------------------
train_imputed$pred_prob <- apply_map(train_imputed$score); test_final$pred_prob <- apply_map(test_final$score)
set.seed(20260813)                                                          # out-of-fold: the map is re-estimated inside each fold
folds <- createFolds(train_imputed$disposition_bin, k = 10)
train_imputed$pred_prob_oof <- NA_real_
for (idx in folds) {
  m <- train_imputed[-idx, ] %>% mutate(g = pmin(pmax(score, POOL_LOW), POOL_HIGH)) %>% group_by(g) %>% summarise(p = mean(disposition_bin), .groups = "drop")
  train_imputed$pred_prob_oof[idx] <- m$p[match(pmin(pmax(train_imputed$score[idx], POOL_LOW), POOL_HIGH), m$g)]
}
cal_tab <- bind_rows(calibration_metrics(train_imputed$disposition_bin, train_imputed$pred_prob, "Training (apparent)"),
                     calibration_metrics(train_imputed$disposition_bin, train_imputed$pred_prob_oof, "Training (out-of-fold)"),
                     calibration_metrics(test_final$disposition_bin, test_final$pred_prob, "Held-out test"))
write.csv(cal_tab, out("calibration_metrics.csv"), row.names = FALSE)

deciles <- function(y, p, label, g = 10) {
  grp <- cut(p, breaks = unique(quantile(p, probs = seq(0, 1, length.out = g + 1), na.rm = TRUE)), include.lowest = TRUE)
  data.frame(y = y, p = p, grp = grp) %>% group_by(grp) %>%
    summarise(n = n(), mean_pred = mean(p), obs = mean(y), obs_low = wilson_ci(sum(y), n())[1], obs_high = wilson_ci(sum(y), n())[2], .groups = "drop") %>% mutate(Partition = label)
}
by_score <- function(score, y, p, label) data.frame(score = score, y = y, p = p) %>% group_by(score) %>%
  summarise(n = n(), predicted = mean(p), observed = mean(y), obs_low = wilson_ci(sum(y), n())[1], obs_high = wilson_ci(sum(y), n())[2], .groups = "drop") %>% mutate(Partition = label)
write.csv(bind_rows(deciles(train_imputed$disposition_bin, train_imputed$pred_prob_oof, "Training (out-of-fold)"),
                    deciles(test_final$disposition_bin, test_final$pred_prob, "Held-out test")), out("supp_fig3_calibration_deciles.csv"), row.names = FALSE)
write.csv(bind_rows(by_score(train_imputed$score, train_imputed$disposition_bin, train_imputed$pred_prob_oof, "Training (out-of-fold)"),
                    by_score(test_final$score, test_final$disposition_bin, test_final$pred_prob, "Held-out test")), out("supp_fig3_calibration_by_score.csv"), row.names = FALSE)

# ---- 4. Performance of the frozen risk bands ------------------------------------
band_of <- function(s) factor(case_when(s <= low_cut ~ "Low", s < high_cut ~ "Intermediate", TRUE ~ "High"), levels = c("Low", "Intermediate", "High"))
band_cal <- bind_rows(train_imputed %>% transmute(Partition = "Training", band = band_of(score), y = disposition_bin, p = pred_prob),
                      test_final %>% transmute(Partition = "Test", band = band_of(score), y = disposition_bin, p = pred_prob)) %>%
  group_by(Partition, band) %>% summarise(N = n(), Observed = sum(y), Observed_rate = mean(y), Expected = sum(p), Mean_predicted = mean(p), OE = sum(y) / sum(p), .groups = "drop_last") %>%
  mutate(Pct_of_cohort = N / sum(N)) %>% ungroup()
write.csv(band_cal, out("risk_band_calibration.csv"), row.names = FALSE)
cat(sprintf("Validation partition: NPV of the low band %.4f | PPV of the high band %.4f\n",
            1 - mean(test_final$disposition_bin[test_final$score <= low_cut]), mean(test_final$disposition_bin[test_final$score >= high_cut])))

# ---- 5. Decision curve: score vs ESI alone vs Admit all / Admit none -------------
esi_fit <- glm(disposition_bin ~ factor(esi), data = train_imputed, family = binomial)     # ESI alone: logistic model fitted on the derivation partition
test_final$prob_esi <- predict(esi_fit, newdata = test_final, type = "response")
ths <- seq(0.01, 0.99, by = 0.01); prev <- mean(test_final$disposition_bin)
nb <- data.frame(threshold = ths, Treat.All = prev - (1 - prev) * ths / (1 - ths), Treat.None = 0,
                 PAnACEA.H.score = sapply(ths, function(t) net_benefit(test_final$disposition_bin, test_final$pred_prob, t)),
                 ESI.alone = sapply(ths, function(t) net_benefit(test_final$disposition_bin, test_final$prob_esi, t))) %>%
  mutate(best_default = pmax(Treat.All, Treat.None), margin_vs_default = PAnACEA.H.score - best_default, margin_vs_esi = PAnACEA.H.score - ESI.alone)
write.csv(nb, out("supp_fig4_decision_curve_net_benefit.csv"), row.names = FALSE)
cat(sprintf("Net benefit exceeds both default strategies from %.2f to %.2f\n", min(nb$threshold[nb$margin_vs_default > 0]), max(nb$threshold[nb$margin_vs_default > 0])))
saveRDS(list(emp_map = emp_map, low_cut = low_cut, high_cut = high_cut, alpha = alpha, beta = beta), file.path(CACHE_DIR, "mapping.rds"))
cat("Done: outputs in", OUT_DIR, "\n")
