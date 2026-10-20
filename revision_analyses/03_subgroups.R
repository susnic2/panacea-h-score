# ==============================================================================
# 03. Performance, equity and case-mix by subgroup
#
# Reproduces (manuscript numbering):
#   Supplementary Table 22  discrimination and calibration by age band, sex, race, ethnicity, insurance status, ED site
#   Supplementary Table 23  magnitude of the disparity within each domain
#   Supplementary Table 24  case-mix by subgroup (TRIPOD+AI item 20b)
#   Supplementary Table 25  distribution of patients across the risk bands and predictive value of each band, by subgroup
#   Supplementary Figure 5  discrimination by subgroup (forest plot)
#   Supplementary Figure 6  calibration by age band
#
# All on the held-out partition. Subgroups with fewer than 500 visits are reported but flagged and not interpreted.
# Run script 01 first. From the repository root:  Rscript revision_analyses/03_subgroups.R
# ==============================================================================
source("revision_analyses/R/setup.R")
map <- readRDS(file.path(CACHE_DIR, "mapping.rds"))
POOL_LOW <- 4; POOL_HIGH <- 33
apply_map <- function(score) map$emp_map$p_map[match(pmin(pmax(score, POOL_LOW), POOL_HIGH), map$emp_map$score_grp)]
test_final$pred_prob <- apply_map(test_final$score)
low_cut <- map$low_cut; high_cut <- map$high_cut
y <- test_final$disposition_bin; p <- test_final$pred_prob; s <- test_final$score
MIN_N <- 500

src <- df[match(test_final$id, df$id), ]                       # subgroup variables come from the full source file
stopifnot(nrow(src) == nrow(test_final), !anyNA(src$id))
subgroups <- list(`Age band` = cut(as.numeric(src$age), c(-Inf, 40, 65, 80, Inf), labels = c("<40", "40-64", "65-79", "80+"), right = FALSE),
                  Sex = src$gender, Race = src$race, Ethnicity = src$ethnicity, Insurance = src$insurance_status, `ED site` = src$dep_name)

# ---- Supplementary Tables 22 and 23 ----------------------------------------------
subgroup_metrics <- function(g, gname) {
  bind_rows(lapply(levels(factor(g)), function(l) {
    i <- which(!is.na(g) & g == l); if (length(i) < 30 || length(unique(y[i])) < 2) return(NULL)
    yi <- y[i]; pi <- p[i]; r <- roc(yi, pi, quiet = TRUE); ci <- ci.auc(r, method = "delong")
    pc <- pmin(pmax(pi, 1e-6), 1 - 1e-6); lp <- qlogis(pc)
    sl <- tryCatch(unname(coef(glm(yi ~ lp, family = binomial))[2]), error = function(e) NA)
    ic <- tryCatch(unname(coef(glm(yi ~ offset(lp), family = binomial))[1]), error = function(e) NA)
    tt <- mean(yi); sel <- pc >= tt                              # net benefit at a threshold equal to the subgroup prevalence
    nb <- sum(yi[sel] == 1) / length(yi) - (sum(yi[sel] == 0) / length(yi)) * (tt / (1 - tt))
    data.frame(Domain = gname, Subgroup = l, N = length(i), Prevalence = mean(yi), AUC = as.numeric(auc(r)), CI_low = ci[1], CI_high = ci[3],
               Cal_intercept = ic, Cal_slope = sl, OE = sum(yi) / sum(pc), Brier = mean((pc - yi)^2), NB_at_prevalence = nb,
               Pct_low_band = mean(s[i] <= low_cut), Pct_high_band = mean(s[i] >= high_cut), Unstable = length(i) < MIN_N)
  }))
}
sub_tab <- bind_rows(lapply(names(subgroups), function(nm) subgroup_metrics(subgroups[[nm]], nm)))
write.csv(sub_tab, out("supp_table22_subgroup_performance.csv"), row.names = FALSE)
disparity <- sub_tab %>% filter(!Unstable) %>% group_by(Domain) %>%
  summarise(AUC_min = min(AUC), AUC_max = max(AUC), AUC_range = max(AUC) - min(AUC), OE_min = min(OE), OE_max = max(OE),
            Intercept_min = min(Cal_intercept), Intercept_max = max(Cal_intercept), .groups = "drop")
write.csv(disparity, out("supp_table23_disparity_summary.csv"), row.names = FALSE)

# ---- Supplementary Table 24: case-mix ---------------------------------------------
pos_rate <- function(x) mean(as.numeric(as.character(x)) > 0, na.rm = TRUE)
case_mix <- bind_rows(lapply(names(subgroups), function(dom) bind_rows(lapply(levels(factor(subgroups[[dom]])), function(l) {
  i <- which(!is.na(subgroups[[dom]]) & subgroups[[dom]] == l); if (length(i) < 30) return(NULL); d <- src[i, ]
  data.frame(Domain = dom, Subgroup = l, N = nrow(d), Admission_rate = mean(d$disposition == "Admit", na.rm = TRUE),
             Median_age = median(as.numeric(d$age), na.rm = TRUE), Pct_ESI_1_2 = mean(as.numeric(as.character(d$esi)) %in% c(1, 2), na.rm = TRUE),
             Median_HR = median(as.numeric(d$triage_vital_hr), na.rm = TRUE), Pct_cardiovascular_meds = pos_rate(d$meds_cardiovascular),
             Pct_analgesic_meds = pos_rate(d$meds_analgesics), Pct_psychiatric_meds = pos_rate(d$meds_psychotherapeuticdrugs),
             Median_n_admissions = median(as.numeric(d$n_admissions), na.rm = TRUE))
}))))
write.csv(case_mix, out("supp_table24_case_mix_by_subgroup.csv"), row.names = FALSE)

# ---- Supplementary Table 25: risk bands by subgroup --------------------------------
band <- factor(ifelse(s <= low_cut, "Low", ifelse(s < high_cut, "Intermediate", "High")), levels = c("Low", "Intermediate", "High"))
one <- function(idx, domain, level) {
  out_row <- data.frame(Domain = domain, Subgroup = level, N = length(idx), Unstable = length(idx) < MIN_N)
  for (b in levels(band)) { ib <- idx[band[idx] == b]
    out_row[[paste0(b, "_N")]] <- length(ib); out_row[[paste0(b, "_pct")]] <- length(ib) / length(idx)
    out_row[[paste0(b, "_admit")]] <- if (length(ib)) mean(y[ib]) else NA }                  # admission rate in the band (Low: 1 - NPV; High: PPV)
  out_row
}
bs <- bind_rows(one(seq_along(y), "All", "All"), bind_rows(lapply(names(subgroups), function(nm) { g <- subgroups[[nm]]
  bind_rows(lapply(levels(factor(g)), function(l) one(which(!is.na(g) & g == l), nm, l))) })))
bs$Low_NPV <- 1 - bs$Low_admit; bs$High_PPV <- bs$High_admit
write.csv(bs, out("supp_table25_risk_bands_by_subgroup.csv"), row.names = FALSE)

# ---- Supplementary Figures 5 and 6 -------------------------------------------------
plot_df <- sub_tab %>% filter(!Unstable) %>% mutate(lab = paste0(Domain, ": ", Subgroup))
p_forest <- ggplot(plot_df, aes(x = AUC, y = reorder(lab, AUC))) +
  geom_vline(xintercept = fast_auc(y, p), linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.25, colour = "steelblue") + geom_point(aes(size = N), colour = "steelblue") +
  scale_size_continuous(range = c(1, 4), guide = "none") + facet_grid(Domain ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(x = "AUC (95% CI)", y = NULL, title = "Discrimination of the PAnACEA-H score by subgroup", subtitle = "Held-out validation partition; dashed line = overall AUC") +
  theme_minimal(base_size = 11) + theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0))
ggsave(out("supp_fig5_subgroup_forest.jpeg"), p_forest, dpi = 300, bg = "white", width = 9, height = 10)
age_band <- subgroups[["Age band"]]
age_cal <- bind_rows(lapply(levels(age_band), function(l) { i <- which(!is.na(age_band) & age_band == l); data.frame(Age = l, score = s[i], y = y[i], p = p[i]) })) %>%
  group_by(Age, score) %>% summarise(n = n(), observed = mean(y), predicted = mean(p), .groups = "drop")
write.csv(age_cal, out("supp_fig6_calibration_by_age.csv"), row.names = FALSE)
p_age <- ggplot(age_cal, aes(x = predicted, y = observed)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
  geom_line(colour = "grey30") + geom_point(aes(size = n), alpha = 0.6, colour = "grey30") + scale_size_continuous(range = c(0.5, 3), guide = "none") +
  facet_wrap(~ Age, nrow = 1) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(x = "Predicted probability", y = "Observed proportion", title = "Calibration by age band", subtitle = "One point per attainable score value; the score awards one point per decade of age") +
  theme_minimal(base_size = 11)
ggsave(out("supp_fig6_calibration_by_age.jpeg"), p_age, dpi = 300, bg = "white", width = 12, height = 3.8)
cat("Done: outputs in", OUT_DIR, "\n")
