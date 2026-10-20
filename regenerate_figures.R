# ==============================================================================
# PAnACEA-H: regenerate the main-text Figures 1-3 and Supplementary Figure 3.
#
# Run from the repository root:   source("regenerate_figures.R")
# Output: ./figures/
#
# Figure 1  Standardized mean differences (SMD, admitted vs discharged) for
#           past medical history, chief complaint and home medications.
#           NOTE: the SMD values below are transcribed from Supplementary
#           Tables 9-11 of the manuscript; they are NOT recomputed from the
#           raw data by this script.
# Figure 2  Top-15 predictors by mean information gain (XGBoost, 5-fold CV).
#           Input: data_for_figures/supplementary_table4_feature_ranking.csv
#           (identical to the output of panacea_h_reproducibility.R).
# Figure 3  ROC curve of the PAnACEA-H score on the held-out test set.
#           Input: reproducibility_output/table3_threshold_performance.csv
#           (sensitivity and specificity at every integer score threshold,
#           which are the exact vertices of the empirical ROC curve of an
#           integer-valued score), and the DeLong AUC and 95% CI stored in
#           reproducibility_output/panacea_h_results.RData.
# Supp. Fig. 3  Calibration of the PAnACEA-H score (held-out test set and
#           out-of-fold derivation partition).
#           Input: data_for_figures/R2_calibration_*.csv, copies of
#           revision_analyses/output/supp_fig3_calibration_by_score.csv and
#           supp_fig3_calibration_deciles.csv, which are written by
#           revision_analyses/01_mapping_calibration_bands.R.
# ==============================================================================

library(dplyr)
library(ggplot2)

OUT <- "figures"
dir.create(OUT, showWarnings = FALSE)

mm_to_in <- function(mm) mm / 25.4

# ==============================================================================
# FIGURE 1 - three-panel SMD bar plot
# ==============================================================================
pmh <- data.frame(
  Category = c("Cardiovascular disease","Endocrine or metabolic disease",
               "Gastrointestinal disease","Kidney disease","Neoplasms",
               "Musculoskeletal disease","Hematological disease",
               "Neurological disease","ENT/Ophthalmological disease","Other",
               "Respiratory disease","Infectious disease","Genital disease",
               "Psychiatric disease","Skin disease","Obstetric disease",
               "Trauma/poisoning"),
  SMD = c(0.67,0.39,0.36,0.36,0.38,0.34,0.23,0.25,0.22,0.19,0.16,0.15,0.11,
          0.11,0.11,0.07,0.05),
  Panel = "A: Past medical history"
)

cc <- data.frame(
  Category = c("Cardiovascular disease","Musculoskeletal disease",
               "Respiratory disease","Neurological disease",
               "ENT/Ophthalmological disease","Trauma/poisoning",
               "Gastrointestinal disease","Genital disease",
               "Endocrine or metabolic disease","Obstetric disease",
               "Kidney disease","Infectious disease","Psychiatric disease",
               "Hematological disease","Other"),
  SMD = c(0.25,0.35,0.23,0.20,0.22,0.19,0.14,0.13,0.08,0.07,0.06,0.04,0.04,
          0.03,0.00),
  Panel = "B: Chief complaint"
)

meds <- data.frame(
  Category = c("Cardiovascular medications","Other","Gastrointestinal medications",
               "Psychiatric medications","Analgesics","Endocrine medications",
               "Neurological medications","Rheumatological medications",
               "Respiratory medications","Antimicrobials","ENT medications",
               "Chemotherapeutics","Hematological medications","Antidotes",
               "Obstetric medications"),
  SMD = c(1.10,0.98,0.91,0.72,0.72,0.76,0.57,0.54,0.56,0.47,0.46,0.19,0.16,
          0.07,0.03),
  Panel = "C: Home medications"
)

fig1_data <- bind_rows(pmh, cc, meds)

# Several category names (e.g. "Cardiovascular disease") repeat across panels
# with different SMDs; a plain shared factor would force one global order on
# all facets. Give each panel its own sort key instead.
fig1_data$sort_key <- paste(fig1_data$Panel, fig1_data$Category, sep = "___")
fig1_data$sort_key <- factor(
  fig1_data$sort_key,
  levels = fig1_data$sort_key[order(fig1_data$Panel, fig1_data$SMD)]
)

p_fig1 <- ggplot(fig1_data, aes(x = SMD, y = sort_key)) +
  geom_col(fill = "steelblue", width = 0.7) +
  facet_wrap(~ Panel, scales = "free_y", ncol = 1) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  labs(x = "Standardized mean difference (admission vs. discharge)",
       y = NULL) +
  theme_minimal(base_size = 10) +
  theme(strip.text = element_text(face = "bold", hjust = 0),
        panel.grid.minor = element_blank())

ggsave(file.path(OUT, "Figure_1.jpeg"), p_fig1, dpi = 300, bg = "white",
       width = mm_to_in(180), height = mm_to_in(220))
cat("Figure 1 saved.\n")

# ==============================================================================
# FIGURE 2 - top-15 information gain (85 mm and 180 mm widths)
# ==============================================================================
feat <- read.csv("data_for_figures/supplementary_table4_feature_ranking.csv")
feat <- head(feat[order(-feat$Mean_Gain), ], 15)

p_fig2 <- ggplot(feat, aes(x = reorder(Feature, Mean_Gain), y = Mean_Gain)) +
  geom_point(size = 2, colour = "steelblue") +
  geom_errorbar(aes(ymin = Mean_Gain - SD_Gain, ymax = Mean_Gain + SD_Gain),
                width = 0.2, linewidth = 0.3) +
  coord_flip() +
  labs(x = NULL, y = "Mean information gain (5-fold CV, +/- SD)") +
  theme_minimal(base_size = 8)

ggsave(file.path(OUT, "Figure_2_singlecolumn_85mm.jpeg"), p_fig2, dpi = 300,
       bg = "white", width = mm_to_in(85), height = mm_to_in(100))
ggsave(file.path(OUT, "Figure_2_doublecolumn_180mm.jpeg"), p_fig2, dpi = 300,
       bg = "white", width = mm_to_in(180), height = mm_to_in(110))
cat("Figure 2 saved (85 mm and 180 mm versions).\n")

# ==============================================================================
# FIGURE 3 - ROC curve, held-out test set (85 mm single-column width)
# ==============================================================================
thr <- read.csv("reproducibility_output/table3_threshold_performance.csv")
res <- new.env(); load("reproducibility_output/panacea_h_results.RData", envir = res)
auc_ci <- as.numeric(res$ci_test)                    # DeLong: lower, AUC, upper

roc_pts <- rbind(data.frame(fpr = 0, sens = 0),      # threshold above the maximum score
                 data.frame(fpr = 1 - thr$Specificity, sens = thr$Sensitivity))
roc_pts <- roc_pts[order(roc_pts$fpr, roc_pts$sens), ]

p_fig3 <- ggplot(roc_pts, aes(fpr, sens)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_line(colour = "steelblue4", linewidth = 0.9) +
  annotate("text", x = 0.97, y = 0.08, hjust = 1, size = 2.8,
           label = sprintf("AUC %.3f\n(95%% CI %.3f-%.3f)", auc_ci[2], auc_ci[1], auc_ci[3])) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank())

ggsave(file.path(OUT, "Figure_3_ROC.jpeg"), p_fig3, dpi = 300, bg = "white",
       width = mm_to_in(85), height = mm_to_in(85))
cat("Figure 3 (ROC) saved.\n")

# ==============================================================================
# SUPPLEMENTARY FIGURE 3 - calibration plot
# ==============================================================================
cal_curve <- read.csv("data_for_figures/R2_calibration_curve_by_score.csv")
cal_dec   <- read.csv("data_for_figures/R2_calibration_plot_data.csv")
names(cal_dec)[names(cal_dec) == "mean_pred"] <- "predicted"

part_levels <- c("Training (out-of-fold)", "Held-out test")
cal_curve$Partition <- factor(cal_curve$Partition, levels = part_levels)
cal_dec$Partition   <- factor(cal_dec$Partition,   levels = part_levels)

p_fig3 <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(data = cal_curve, aes(x = predicted, y = observed),
            colour = "grey30", linewidth = 0.5) +
  geom_point(data = cal_curve, aes(x = predicted, y = observed, size = n),
             colour = "grey30", alpha = 0.5) +
  scale_size_continuous(range = c(0.4, 3), guide = "none") +
  geom_errorbar(data = cal_dec,
                aes(x = predicted, ymin = obs_low, ymax = obs_high),
                width = 0.01, colour = "steelblue") +
  geom_point(data = cal_dec, aes(x = predicted, y = obs),
             size = 1.6, colour = "steelblue") +
  facet_wrap(~ Partition) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability of admission",
       y = "Observed proportion admitted",
       title = "Calibration of the PAnACEA-H score") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 11))

ggsave(file.path(OUT, "Supplementary_Figure_3_calibration.jpeg"), p_fig3, dpi = 300,
       bg = "white", width = mm_to_in(180), height = mm_to_in(100))
cat("Supplementary Figure 3 (calibration) saved.\n")

cat("\nAll figures written to:", normalizePath(OUT), "\n")
