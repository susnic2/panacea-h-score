# ==============================================================================
# 02. Comparator models, formal comparisons and net benefit
#
# Reproduces (manuscript numbering):
#   Supplementary Table 19  net benefit of the score and of each comparator
#   Supplementary Table 20  AUC (95% CI) and calibration of every model
#   Supplementary Table 21  formal comparisons against the score (DeLong; bootstrap of the difference)
#   START figures quoted in the text (AUC of the raw score, observed-to-expected ratio)
#
# Models, all evaluated once on the held-out partition:
#   PAnACEA-H     probabilities from the score-to-probability map (script 01)
#   ESI alone     logistic model on the five ESI levels, fitted on the derivation partition
#   START         published integer weights (Dinh et al., BMC Emerg Med 2016) and published probability
#                 bands; and the same variable set refitted by logistic regression on the derivation partition
#   Hong et al.   low-dimensional model rebuilt from its description: 14 patient and hospital-use variables
#                 plus the 48 medication-class counts (62 inputs), XGBoost depth 6, eta 0.3, 100 rounds
#   XGBoost       the model of the Methods: full 545-variable pool, depth 3, eta 0.1, 100 rounds, median/mode imputation
#
# All comparisons use the visits for which every comparator can be computed (START needs ESI, which is
# missing for 779 visits). Run script 01 first. From the repository root:
#   Rscript revision_analyses/02_comparators.R
# ==============================================================================
source("revision_analyses/R/setup.R")
map <- readRDS(file.path(CACHE_DIR, "mapping.rds"))
POOL_LOW <- 4; POOL_HIGH <- 33
apply_map <- function(score) map$emp_map$p_map[match(pmin(pmax(score, POOL_LOW), POOL_HIGH), map$emp_map$score_grp)]
test_final$pred_prob <- apply_map(test_final$score)
N_BOOT <- 2000

# ---- 1. ESI alone -------------------------------------------------------------
esi_fit <- glm(disposition_bin ~ factor(esi), data = train_imputed, family = binomial)
test_final$prob_esi <- predict(esi_fit, newdata = test_final, type = "response")

# ---- 2. XGBoost, model of the Methods -------------------------------------------
imp <- impute_median_mode(train_data_raw, test_data_raw)
xtr <- prepare_xgb_data(imp$train); xte <- prepare_xgb_data(imp$test)
cols <- intersect(setdiff(names(xtr), "disposition"), setdiff(names(xte), "disposition"))
stopifnot(length(cols) == 545)
set.seed(2024)
xgb_methods <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3, eta = 0.1, nthread = 2),
                         data = xgb.DMatrix(as.matrix(xtr[, cols]), label = xtr$disposition), nrounds = 100, verbose = 0)
test_final$prob_xgb <- predict(xgb_methods, xgb.DMatrix(as.matrix(xte[, cols]), label = xte$disposition))
cat(sprintf("XGBoost (Methods model) held-out AUC: %.5f (reported 0.91389)\n", fast_auc(y_test, test_final$prob_xgb)))
stopifnot(abs(fast_auc(y_test, test_final$prob_xgb) - 0.91389) < 5e-4)

# ---- 3. Low-dimensional model of Hong et al. (native handling of missing values) ----
hong_named <- c("esi", "age", "gender", "maritalstatus", "employstatus", "insurance_status", "race", "ethnicity",
                "lang", "religion", "n_edvisits", "n_admissions", "previousdispo", "n_surgeries")
hong_vars <- c(hong_named, grep("^meds", names(df), value = TRUE))
stopifnot(length(hong_vars) == 62)
tr_h <- df[match(train_data_raw$id, df$id), hong_vars]; te_h <- df[match(test_data_raw$id, df$id), hong_vars]
encode <- function(a, b) {
  for (v in names(a)) if (is.factor(a[[v]]) || is.character(a[[v]])) {
    lev <- levels(factor(a[[v]])); a[[v]] <- as.numeric(factor(a[[v]], levels = lev)); b[[v]] <- as.numeric(factor(b[[v]], levels = lev))
  } else { a[[v]] <- as.numeric(a[[v]]); b[[v]] <- as.numeric(b[[v]]) }
  list(a = as.matrix(a), b = as.matrix(b))
}
enc <- encode(tr_h, te_h)
set.seed(2024)
xgb_hong <- xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", max_depth = 6, eta = 0.3, nthread = 2),
                      data = xgb.DMatrix(enc$a, label = y_train), nrounds = 100, verbose = 0)
test_final$prob_hong <- predict(xgb_hong, xgb.DMatrix(enc$b))
cat(sprintf("Hong low-dimensional model held-out AUC: %.4f (source study: 0.910)\n", fast_auc(y_test, test_final$prob_hong)))

# ---- 4. START (published weights and refitted variable set) ---------------------
# Limits of transportability, to be stated wherever START is reported: START uses the Australasian Triage Scale
# (ESI is used as a surrogate); it assumes one presenting problem per visit, whereas the source dataset has 200
# non-exclusive chief-complaint flags (a priority rule assigns one); its "admission in the previous 30 days"
# is proxied by "at least one admission in the previous year"; and its outcome includes short-stay and transfers.
start_age_points <- function(age) case_when(age < 20 ~ 0, age < 40 ~ 1, age < 60 ~ 3, age < 80 ~ 6, TRUE ~ 9)
start_triage_points <- c("1" = 24, "2" = 16, "3" = 11, "4" = 5, "5" = 0)
start_ambulance_points <- function(mode) ifelse(grepl("ambulance", as.character(mode), ignore.case = TRUE), 4, 0)
start_prevadm_points <- function(n_adm) ifelse(n_adm >= 1, 3, 0)
start_hour_points <- function(bin) ifelse(as.character(bin) %in% c("07-10", "11-14", "15-18"), 1, 0)
start_pp_points <- c(abdominal_gi = 2, cardiovascular = -3, general = 0, febrile = 3, injury = -4, respiratory = 0, musculoskeletal = -3,
                     neurological = -1, mental_health = -2, toxicological = -2, ent_eye_head_neck = -6, administrative = -3,
                     genitourinary = -1, social = 1, endocrine = 0, obs_gyn = -3, skin_allergy = -2, other_medical = 5)
# chief-complaint flags -> START presenting-problem categories
cc_map <- list(
  abdominal_gi = c("cc_abdominalcramping","cc_abdominaldistention","cc_abdominalpain",
                   "cc_constipation","cc_diarrhea","cc_emesis","cc_epigastricpain","cc_gibleeding",
                   "cc_giproblem","cc_nausea","cc_rectalbleeding","cc_rectalpain",
                   "cc_swallowedforeignbody","cc_flankpain"),
  cardiovascular = c("cc_cardiacarrest","cc_chestpain","cc_chesttightness","cc_edema",
                     "cc_hypertension","cc_hypotension","cc_irregularheartbeat","cc_palpitations",
                     "cc_rapidheartrate","cc_tachycardia","cc_legswelling"),
  febrile = c("cc_fever","cc_fever-75yearsorolder","cc_fever-9weeksto74years",
              "cc_feverimmunocompromised","cc_chills","cc_influenza","cc_uri",
              "cc_urinarytractinfection","cc_cellulitis","cc_abscess","cc_woundinfection",
              "cc_follow-upcellulitis"),
  injury = c("cc_animalbite","cc_ankleinjury","cc_arminjury","cc_assaultvictim",
             "cc_burn","cc_extremitylaceration","cc_eyeinjury","cc_facialinjury",
             "cc_faciallaceration","cc_fall","cc_fall>65","cc_fingerinjury","cc_footinjury",
             "cc_fulltrauma","cc_handinjury","cc_headinjury","cc_headlaceration",
             "cc_kneeinjury","cc_laceration","cc_leginjury","cc_modifiedtrauma",
             "cc_motorcyclecrash","cc_motorvehiclecrash","cc_multiplefalls","cc_ribinjury",
             "cc_shoulderinjury","cc_thumbinjury","cc_toeinjury","cc_trauma",
             "cc_wristinjury","cc_insectbite","cc_tickremoval"),
  respiratory = c("cc_asthma","cc_breathingdifficulty","cc_breathingproblem",
                  "cc_cough","cc_dyspnea","cc_hemoptysis","cc_respiratorydistress",
                  "cc_shortnessofbreath","cc_wheezing","cc_coldlikesymptoms","cc_nasalcongestion"),
  musculoskeletal = c("cc_anklepain","cc_armpain","cc_armswelling","cc_backpain",
                      "cc_elbowpain","cc_fingerpain","cc_fingerswelling","cc_footpain",
                      "cc_footswelling","cc_generalizedbodyaches","cc_groinpain","cc_handpain",
                      "cc_hippain","cc_jointswelling","cc_kneepain","cc_legpain","cc_neckpain",
                      "cc_ribpain","cc_shoulderpain","cc_toepain","cc_wristpain","cc_jawpain"),
  neurological = c("cc_alteredmentalstatus","cc_confusion","cc_dizziness",
                   "cc_extremityweakness","cc_headache","cc_headache-newonsetornewsymptoms",
                   "cc_headache-recurrentorknowndxmigraines","cc_headachere-evaluation",
                   "cc_lossofconsciousness","cc_migraine","cc_nearsyncope","cc_neurologicproblem",
                   "cc_numbness","cc_seizure-newonset","cc_seizure-priorhxof","cc_seizures",
                   "cc_strokealert","cc_syncope","cc_unresponsive","cc_weakness","cc_lethargy"),
  mental_health = c("cc_addictionproblem","cc_agitation","cc_alcoholproblem",
                    "cc_anxiety","cc_depression","cc_detoxevaluation","cc_drug/alcoholassessment",
                    "cc_drugproblem","cc_hallucinations","cc_homicidal","cc_panicattack",
                    "cc_psychiatricevaluation","cc_psychoticsymptoms","cc_suicidal",
                    "cc_withdrawal-alcohol"),
  toxicological = c("cc_alcoholintoxication","cc_ingestion","cc_overdose-accidental",
                    "cc_overdose-intentional","cc_poisoning"),
  ent_eye_head_neck = c("cc_blurredvision","cc_conjunctivitis","cc_dentalpain",
                        "cc_earpain","cc_earproblem","cc_epistaxis","cc_eyepain","cc_eyeproblem",
                        "cc_eyeredness","cc_facialpain","cc_facialswelling","cc_foreignbodyineye",
                        "cc_oralswelling","cc_otalgia","cc_sinusproblem","cc_sorethroat"),
  administrative = c("cc_abnormallab","cc_medicalscreening","cc_medicationrefill",
                     "cc_suture/stapleremoval","cc_woundcheck","cc_woundre-evaluation",
                     "cc_post-opproblem","cc_medicationproblem"),
  genitourinary = c("cc_dysuria","cc_hematuria","cc_maleguproblem",
                    "cc_testiclepain","cc_urinaryfrequency","cc_urinaryretention"),
  social = character(0),
  endocrine = c("cc_decreasedbloodsugar-symptomatic","cc_dehydration",
                "cc_elevatedbloodsugar-nosymptoms","cc_elevatedbloodsugar-symptomatic",
                "cc_hyperglycemia"),
  obs_gyn = c("cc_abdominalpainpregnant","cc_breastpain","cc_femaleguproblem",
              "cc_pelvicpain","cc_vaginalbleeding","cc_vaginaldischarge","cc_vaginalpain"),
  skin_allergy = c("cc_allergicreaction","cc_bleeding/bruising","cc_cyst","cc_mass",
                   "cc_rash","cc_skinirritation","cc_skinproblem"),
  other_medical = c("cc_medicalproblem","cc_sicklecellpain","cc_bodyfluidexposure",
                    "cc_exposuretostd","cc_stdcheck"),
  general = c("cc_fatigue","cc_pain","cc_other")
)

START_PRIORITY <- c("cardiovascular", "respiratory", "neurological", "toxicological", "febrile", "abdominal_gi", "injury", "mental_health",
                    "endocrine", "obs_gyn", "genitourinary", "musculoskeletal", "ent_eye_head_neck", "skin_allergy", "other_medical", "social",
                    "administrative", "general")                 # clinical priority, fixed in advance, used to assign ONE category per visit
start_presenting_problem <- function(data) {
  flag <- matrix(FALSE, nrow(data), length(cc_map), dimnames = list(NULL, names(cc_map)))
  for (k in names(cc_map)) {
    cols <- intersect(cc_map[[k]], names(data))
    if (length(cols)) { m <- sapply(data[, cols, drop = FALSE], function(x) as.numeric(as.character(x)) > 0); if (is.null(dim(m))) m <- matrix(m, ncol = 1); flag[, k] <- rowSums(m, na.rm = TRUE) > 0 }
  }
  primary <- rep(NA_character_, nrow(data))
  for (k in START_PRIORITY) { i <- is.na(primary) & flag[, k]; primary[i] <- k }
  primary[is.na(primary)] <- "general"
  list(category = primary, points = unname(start_pp_points[primary]))
}
start_src_te <- df[match(test_final$id, df$id), ]; start_src_tr <- df[match(train_imputed$id, df$id), ]
pp_te <- start_presenting_problem(start_src_te); pp_tr <- start_presenting_problem(start_src_tr)
test_final$start_score <- start_age_points(as.numeric(start_src_te$age)) + unname(start_triage_points[as.character(start_src_te$esi)]) +
  start_ambulance_points(start_src_te$arrivalmode) + start_prevadm_points(as.numeric(start_src_te$n_admissions)) +
  start_hour_points(start_src_te$arrivalhour_bin) + pp_te$points                         # NA where ESI is missing
# published probability bands; a missing score stays missing (it must NOT fall into the top band)
start_published_prob <- function(s) { p <- case_when(s < 1 ~ 0.03, s <= 10 ~ 0.14, s <= 20 ~ 0.47, s <= 30 ~ 0.81, s <= 40 ~ 0.96, s > 40 ~ 0.99, TRUE ~ NA_real_); p[is.na(s)] <- NA_real_; p }
test_final$prob_start <- start_published_prob(test_final$start_score)
fit_df <- data.frame(y = y_train, age_pts = start_age_points(as.numeric(start_src_tr$age)), esi = factor(start_src_tr$esi),
                     amb = start_ambulance_points(start_src_tr$arrivalmode), prev = start_prevadm_points(as.numeric(start_src_tr$n_admissions)),
                     hour = start_hour_points(start_src_tr$arrivalhour_bin), pp = factor(pp_tr$category))
new_df <- data.frame(age_pts = start_age_points(as.numeric(start_src_te$age)), esi = factor(start_src_te$esi, levels = levels(fit_df$esi)),
                     amb = start_ambulance_points(start_src_te$arrivalmode), prev = start_prevadm_points(as.numeric(start_src_te$n_admissions)),
                     hour = start_hour_points(start_src_te$arrivalhour_bin), pp = factor(pp_te$category, levels = levels(fit_df$pp)))
start_refit <- glm(y ~ age_pts + esi + amb + prev + hour + pp, data = fit_df, family = binomial)
test_final$prob_startrefit <- predict(start_refit, newdata = new_df, type = "response")

# ---- 5. All models on the same visits ---------------------------------------------
models <- list("PAnACEA-H score" = test_final$pred_prob, "XGBoost (Methods model)" = test_final$prob_xgb, "Hong low-dimensional" = test_final$prob_hong,
               "START (published weights)" = test_final$prob_start, "START (refitted here)" = test_final$prob_startrefit, "ESI alone" = test_final$prob_esi)
complete_idx <- Reduce(`&`, lapply(models, function(p) !is.na(p)))
cat(sprintf("Visits with every comparator computable: %d of %d (%.2f%%)\n", sum(complete_idx), length(complete_idx), 100 * mean(complete_idx)))
stopifnot(sum(complete_idx) == 167366)
y_c <- y_test[complete_idx]; models_c <- lapply(models, function(p) p[complete_idx])
cat(sprintf("START, raw score: AUC %.4f on these visits (text: 0.836)\n", fast_auc(y_c, test_final$start_score[complete_idx])))

auc_tab <- bind_rows(lapply(names(models_c), function(nm) { r <- roc(y_c, models_c[[nm]], quiet = TRUE); ci <- ci.auc(r, method = "delong")
  data.frame(Model = nm, N = length(y_c), AUC = as.numeric(auc(r)), CI_low = ci[1], CI_high = ci[3]) }))
cal_tab <- bind_rows(lapply(names(models_c), function(nm) cbind(Model = nm, calibration_metrics(y_c, models_c[[nm]], "Held-out test, complete cases"))))
write.csv(auc_tab, out("supp_table20_auc_all_models.csv"), row.names = FALSE); write.csv(cal_tab, out("supp_table20_calibration_all_models.csv"), row.names = FALSE)

# formal comparisons against the score: DeLong, and a bootstrap of the AUC difference
rocs <- lapply(models_c, function(p) roc(y_c, p, quiet = TRUE)); ref <- rocs[["PAnACEA-H score"]]
comparators <- setdiff(names(models_c), "PAnACEA-H score")
delong_tab <- bind_rows(lapply(comparators, function(nm) { tt <- roc.test(ref, rocs[[nm]], method = "delong")
  data.frame(Comparator = nm, AUC_score = as.numeric(auc(ref)), AUC_comparator = as.numeric(auc(rocs[[nm]])),
             Delta = as.numeric(auc(ref)) - as.numeric(auc(rocs[[nm]])), Z = unname(tt$statistic), p = tt$p.value) }))
set.seed(20260813); n <- length(y_c)
boot <- matrix(NA_real_, N_BOOT, length(comparators), dimnames = list(NULL, comparators))
for (b in seq_len(N_BOOT)) {
  i <- sample.int(n, n, replace = TRUE); a_ref <- fast_auc(y_c[i], models_c[["PAnACEA-H score"]][i])
  for (nm in comparators) boot[b, nm] <- a_ref - fast_auc(y_c[i], models_c[[nm]][i])
}
delong_tab$Boot_CI_low <- apply(boot, 2, quantile, 0.025)[delong_tab$Comparator]; delong_tab$Boot_CI_high <- apply(boot, 2, quantile, 0.975)[delong_tab$Comparator]
write.csv(delong_tab, out("supp_table21_formal_comparisons.csv"), row.names = FALSE)

# ---- 6. Net benefit of every model (Supplementary Table 19) -----------------------
ths <- seq(0.01, 0.99, by = 0.01); prev <- mean(y_c)
nb <- data.frame(threshold = ths, `Admit all` = prev - (1 - prev) * ths / (1 - ths), `Admit none` = 0, check.names = FALSE)
for (nm in names(models_c)) nb[[nm]] <- sapply(ths, function(t) net_benefit(y_c, models_c[[nm]], t))
nb$n <- length(y_c)
write.csv(nb, out("supp_table19_net_benefit_all_models.csv"), row.names = FALSE)
print(nb[nb$threshold %in% c(0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.70, 0.90), ], digits = 4, row.names = FALSE)
cat("Done: outputs in", OUT_DIR, "\n")
