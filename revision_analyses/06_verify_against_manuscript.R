# ==============================================================================
# 06. Check the output of scripts 01-05 against the Supplementary Tables as printed
#
# expected/supp_table_NN.csv holds the Supplementary Tables exactly as they appear in
# the Supplementary Material (formatted strings). For every cell this script rebuilds
# the printed table from the CSV files in revision_analyses/output/ and compares the
# two at the precision shown: a printed 0.8799 must equal the output rounded to four
# decimals, a printed 59.7% must equal 100 x the output rounded to one decimal.
#
# From the repository root:   Rscript revision_analyses/06_verify_against_manuscript.R
# Exit status 0 = every cell agrees; 1 = at least one cell differs (listed below).
#
# Supplementary Tables 2, 3, 9-12, 15 and 16 are produced by panacea_h_reproducibility.R
# and are not part of this folder.
# ==============================================================================
OUT_DIR <- "revision_analyses/output"; EXP_DIR <- "revision_analyses/expected"
rd  <- function(f) read.csv(file.path(OUT_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)
exp_ <- function(f) read.csv(file.path(EXP_DIR, f), check.names = FALSE, colClasses = "character")

# ---- parsing of a printed cell -------------------------------------------------------
parse_printed <- function(s) {
  s <- trimws(gsub("\u2212", "-", s))
  if (s %in% c("", "-", "nan", "NaN", "NA")) return(list(type = "na"))
  if (grepl("^<\\s*0\\.[0-9]+$", s)) return(list(type = "lt", bound = as.numeric(sub("^<\\s*", "", s))))
  pct <- grepl("%$", s); t <- gsub("[,%+]", "", s)
  if (grepl("^-?[0-9]*\\.?[0-9]+[eE]-?[0-9]+$", t)) return(list(type = "sci", value = as.numeric(t)))
  if (grepl("^-?[0-9]*\\.?[0-9]+$", t)) {
    dec <- if (grepl("\\.", t)) nchar(sub("^.*\\.", "", t)) else 0
    return(list(type = "num", value = as.numeric(t), dec = dec, scale = if (pct) 100 else 1))
  }
  list(type = "text", value = s)
}
cell_ok <- function(printed, v) {
  if (!is.na(v) && identical(trimws(printed), trimws(as.character(v)))) return(TRUE)     # identical strings, e.g. "<40" or "5%"
  p <- parse_printed(printed)
  if (p$type == "na")   return(is.na(v) || identical(as.character(v), "NaN"))
  if (p$type == "lt")   return(!is.na(suppressWarnings(as.numeric(v))) && as.numeric(v) < p$bound)
  if (p$type == "text") { if (is.na(v)) return(FALSE)
                          return(toupper(trimws(as.character(v))) == toupper(sub(" \\*$", "", p$value))) }   # trailing " *" is a footnote marker
  x <- suppressWarnings(as.numeric(v)); if (is.na(x)) return(FALSE)
  if (p$type == "sci") return(abs(x - p$value) <= 0.0011 * abs(p$value))                   # three significant digits
  abs(x * p$scale - p$value) <= 0.5 * 10^(-p$dec) + 1e-9
}

results <- list()
check <- function(label, expected_file, view, text_cols = TRUE) {
  e <- exp_(expected_file)
  if (!identical(dim(e), dim(view))) {
    results[[label]] <<- data.frame(table = label, cells = NA, bad = NA, note = sprintf("shape: printed %dx%d, rebuilt %dx%d", nrow(e), ncol(e), nrow(view), ncol(view)))
    return(invisible()) }
  bad <- character(0); n <- 0
  for (j in seq_len(ncol(e))) for (i in seq_len(nrow(e))) {
    pr <- e[i, j]; if (!text_cols && parse_printed(pr)$type == "text") next
    n <- n + 1
    if (!cell_ok(pr, view[i, j])) bad <- c(bad, sprintf("[%d,%s] printed '%s' vs %s", i, names(e)[j], pr, format(view[i, j], digits = 7)))
  }
  results[[label]] <<- data.frame(table = label, cells = n, bad = length(bad), note = paste(head(bad, 3), collapse = " | "))
  if (length(bad)) cat(sprintf("\n%s\n  %s\n", label, paste(head(bad, 8), collapse = "\n  ")))
}

# ---- tables whose printed layout equals the CSV ---------------------------------------
check("Supplementary Table 1",  "supp_table_01.csv", rd("supp_table1_participant_flow.csv"))
check("Supplementary Table 4",  "supp_table_04.csv", rd("supp_table4_predictor_set_sensitivity.csv"), text_cols = FALSE)  # variant labels are translated in the printed table
check("Supplementary Table 5",  "supp_table_05.csv", rd("supp_table5_medication_variables.csv"))
check("Supplementary Table 6",  "supp_table_06.csv", rd("supp_table6_optimism_corrected.csv"))
check("Supplementary Table 7",  "supp_table_07.csv", rd("supp_table7_selection_stability.csv"))
check("Supplementary Table 8",  "supp_table_08.csv", rd("supp_table8_missingness_by_variable.csv"))
check("Supplementary Table 13", "supp_table_13.csv", rd("supp_table13_importance_comparison.csv"))
check("Supplementary Table 17", "supp_table_17.csv", rd("supp_table17_score_to_probability.csv"))
e21 <- exp_("supp_table_21.csv"); t21 <- rd("supp_table21_formal_comparisons.csv")
check("Supplementary Table 21", "supp_table_21.csv", t21[match(e21$Comparator, t21$Comparator), ])
check("Supplementary Table 22", "supp_table_22.csv", rd("supp_table22_subgroup_performance.csv"))
check("Supplementary Table 23", "supp_table_23.csv", rd("supp_table23_disparity_summary.csv"))
check("Supplementary Table 24", "supp_table_24.csv", rd("supp_table24_case_mix_by_subgroup.csv"))

# ---- tables printed in a different layout from the CSV ----------------------------------
t18 <- rd("supp_table18_table3_threshold_metrics_both_partitions.csv")
check("Supplementary Table 18", "supp_table_18.csv", subset(t18, Partition == "Training")[, c("Threshold", "Sensitivity", "Specificity", "PPV", "NPV", "Accuracy")])

e19 <- exp_("supp_table_19.csv"); t19 <- rd("supp_table19_net_benefit_all_models.csv")
v19 <- t19[match(round(as.numeric(sub("%", "", e19$Threshold)) / 100, 2), round(t19$threshold, 2)), ]
v19 <- data.frame(Threshold = e19$Threshold, v19[, setdiff(names(e19), "Threshold"), drop = FALSE], check.names = FALSE)
check("Supplementary Table 19", "supp_table_19.csv", v19)

e20 <- exp_("supp_table_20.csv"); a20 <- rd("supp_table20_auc_all_models.csv"); c20 <- rd("supp_table20_calibration_all_models.csv")
c20 <- c20[grepl("complete", c20$Partition), ]
key20 <- function(lab) { lab <- sub(" \\(Methods model.*$", " (Methods model)", lab); lab }
v20 <- do.call(rbind, lapply(key20(e20$Model), function(m) {
  a <- a20[a20$Model == m, ]; cc <- c20[c20$Model == m, ]
  if (nrow(cc) == 0) cc <- data.frame(OE_ratio = NA, Cal_intercept = NA, Cal_slope = NA, Brier = NA, Brier_scaled = NA)   # raw START score: no probabilities, no calibration
  data.frame(Model = m, AUC = a$AUC, CI_low = a$CI_low, CI_high = a$CI_high, OE_ratio = cc$OE_ratio, Cal_intercept = cc$Cal_intercept,
             Cal_slope = cc$Cal_slope, Brier = cc$Brier, Brier_scaled = cc$Brier_scaled) }))
v20$Model <- e20$Model
check("Supplementary Table 20", "supp_table_20.csv", v20)

e25 <- exp_("supp_table_25.csv"); t25 <- rd("supp_table25_risk_bands_by_subgroup.csv")
v25 <- data.frame(Domain = ifelse(t25$Domain == "All", "", t25$Domain), Subgroup = t25$Subgroup, N = t25$N,
                  `Low: share` = t25$Low_pct, `Low: NPV` = t25$Low_NPV, `Intermediate: share` = t25$Intermediate_pct,
                  `Intermediate: admission rate` = t25$Intermediate_admit, `High: share` = t25$High_pct, `High: PPV` = t25$High_PPV, check.names = FALSE)
check("Supplementary Table 25", "supp_table_25.csv", v25)

e26a <- exp_("supp_table_26a.csv"); a26 <- rd("supp_table26_panelA_medication_missingness_by_arm.csv")
# the CSV stores percentages, the printed table shows them with a % sign
cls <- c(Cardiovascular = "meds_cardiovascular", Gastrointestinal = "meds_gastrointestinal", Analgesic = "meds_analgesics", Psychotherapeutic = "meds_psychotherapeuticdrugs")
v26a <- do.call(rbind, lapply(seq_len(nrow(e26a)), function(i) {
  k <- cls[match(e26a[i, 1], names(cls))]; r <- a26[a26$disposition == e26a[i, 2], ]
  data.frame(Class = e26a[i, 1], Arm = e26a[i, 2], N = r$n, miss = r[[paste0(k, "__pct_missing")]] / 100, nz = r[[paste0(k, "__pct_nonzero_amongNonMissing")]] / 100) }))
check("Supplementary Table 26, panel A", "supp_table_26a.csv", setNames(v26a, names(e26a)))

b26 <- rd("supp_table26_panelB_triage_only_score.csv")
check("Supplementary Table 26, panel B", "supp_table_26b.csv",
      data.frame(b26[, 1], b26$AUC_test, c(NA, b26$AUC_test[-1] - b26$AUC_test[1]), c("reference", ifelse(b26$DeLong_p[-1] < 0.001, "<0.001", b26$DeLong_p[-1]))), text_cols = FALSE)

# Supplementary Table 14 is assembled from three files
e14 <- exp_("supp_table_14.csv"); ag <- rd("supp_table14_age_encoding.csv"); hr <- rd("supp_table14_hr_cutoff_sensitivity.csv"); rw <- rd("supp_table14_hr_weight_reoptimised.csv")
v14 <- do.call(rbind, lapply(seq_len(nrow(e14)), function(i) {
  if (e14[i, 1] == "Age encoding") {
    r <- ag[if (grepl("^Decade", e14[i, 2])) 1 else 2, ]
    data.frame(a = e14[i, 1], b = e14[i, 2], c = NA, d = r$Test_AUC, e = if (grepl("^Decade", e14[i, 2])) "reference" else r$Difference, f = NA, g = NA)
  } else {
    cut <- as.numeric(sub(" bpm.*$", "", e14[i, 2])); r <- hr[abs(hr$HR_cutoff - cut) < 1e-6, ]; w <- rw$weight[abs(rw$cutoff - cut) < 1e-6]
    data.frame(a = e14[i, 1], b = e14[i, 2], c = r$AUC_train, d = r$AUC_test, e = if (grepl("published", e14[i, 2])) "reference" else r$Delta_AUC_test_vs_published,
               f = r$Pct_train_above_cutoff, g = if (length(w)) w else NA) }
}))
check("Supplementary Table 14", "supp_table_14.csv", setNames(v14, names(e14)), text_cols = FALSE)

# ---- summary ----------------------------------------------------------------------------
res <- do.call(rbind, results); rownames(res) <- NULL
cat("\n=== Agreement with the Supplementary Tables as printed ===\n"); print(res[, 1:3], row.names = FALSE)
failed <- sum(res$bad, na.rm = TRUE) + sum(is.na(res$cells))
cat(sprintf("\nCells checked: %d | cells that differ: %d | tables with a layout problem: %d\n",
            sum(res$cells, na.rm = TRUE), sum(res$bad, na.rm = TRUE), sum(is.na(res$cells))))
if (failed > 0) { cat("VERIFICATION FAILED\n"); quit(status = 1) } else cat("VERIFICATION PASSED\n")
