# Analyses added during peer review

Scripts for the comparator, calibration, subgroup, optimism and sensitivity analyses reported in the
revised manuscript and its Supplementary Material. They complement `panacea_h_reproducibility.R`, which
is left unchanged and reproduces the derivation and validation of the score.

Run everything from the **repository root**, with the source dataset (`5v_cleandf.RData`, see the main
README) either in the root or pointed to by `PANACEA_DATA`:

```bash
PANACEA_DATA=/path/to/5v_cleandf.RData bash revision_analyses/run_all.sh
```

or one script at a time, in order (`Rscript revision_analyses/01_mapping_calibration_bands.R`, and so on).
`R/setup.R` is sourced by every script. It rebuilds the 70/30 split exactly as the main pipeline does,
re-derives the score on the derivation partition (multiple imputation with seed 789, Youden cutoffs,
greedy weights), stops if the weights are not the published ones (3, 1, 5, 5, 4, 2, 2) and caches the
objects that take long to compute in `cache/`. Results go to `output/`.

## Scripts and what they reproduce

Numbering follows the revised manuscript.

| Script | Reproduces | Run time* |
|---|---|---|
| `01_mapping_calibration_bands.R` | Supp. Table 17 (score to probability), Supp. Table 18 and Table 3 (threshold metrics), risk-band boundaries and performance, data behind Supp. Fig. 3 (calibration) and Supp. Fig. 4 (decision curve) | 1-5 min |
| `02_comparators.R` | Supp. Tables 19-21: net benefit, AUC and calibration, DeLong and bootstrap comparisons for ESI alone, START (published weights; refitted; raw score), the low-dimensional model of Hong et al. and the XGBoost model of the Methods | 1-5 min |
| `03_subgroups.R` | Supp. Tables 22-25 (subgroup performance, disparity summary, case-mix, risk bands by subgroup), Supp. Figs. 5-6 | 1-5 min |
| `04_optimism_and_selection.R` | Supp. Tables 6-7: Harrell optimism correction (100 bootstrap replicates, two versions) and stability of automated predictor selection | about 6 min on 8 cores (`PANACEA_CORES`) |
| `05_sensitivity_analyses.R` | Supp. Tables 1, 4, 5, 8, 13, 14, 26 and the ranking with and without imputation reported in the Methods | about 10 min |
| `06_verify_against_manuscript.R` | Compares the output with the Supplementary Tables as printed (`expected/`) | seconds |

\*Measured on an 18-core Mac (the whole folder, from an empty `cache/` and `output/`, took 21 minutes); 04 needs several GB of memory per core.

`06_verify_against_manuscript.R` rebuilds each printed table from `output/` and compares every cell at the
precision shown (a printed 0.8799 must equal the output rounded to four decimals; a printed 59.7% must equal
100 times the output rounded to one decimal). It exits with status 0 if all cells agree. Supp. Tables 2, 3,
9-12, 15 and 16 are produced by the main pipeline, not by this folder.

## Choices that matter for interpretation

* **Comparator set.** All comparisons in `02` use the 167,366 held-out visits (99.5%) for which every
  comparator can be computed; START requires the Emergency Severity Index, which is missing for 779 visits.
* **Hong et al. model.** Rebuilt from its published description (14 patient and hospital-use variables plus
  the 48 medication-class counts, XGBoost depth 6, learning rate 0.3, 100 rounds), not from the authors' code.
  It reaches an AUC of 0.914 here against the 0.910 they report.
* **START (published weights)** is evaluated twice: through the published six-level probability bands (the
  version used for calibration and net benefit) and as the raw integer score (the AUC quoted in the text; it has
  its own row in Supp. Tables 20 and 21, without calibration measures because it is not a probability).
* **Optimism correction (`04`).** Missing values are imputed with the median of each bootstrap sample because
  repeating multiple imputation 100 times is not practical; the apparent AUC of this pipeline (0.8808) is
  therefore slightly different from the one on the multiply imputed derivation partition (0.8826). The
  "full pipeline" version replaces the selection rule of the paper (importance, non-redundancy, clinical
  interpretability) by its automated part, the seven variables with the highest information gain among the
  eligible candidates. This rule does not reproduce the published predictor set: gastrointestinal medications
  are selected in every replicate and triage heart rate in none (Supp. Table 7).
* **Gain, permutation importance and SHAP (`05`).** The intercept column of the SHAP contributions is the last
  one and is dropped by position, because the xgboost package names it `BIAS` or `(Intercept)` depending on
  the version.

## Software

R 4.6.1 with dplyr 1.2.1, tidyr 1.3.2, caret 7.0.1, xgboost 3.2.1.1, pROC 1.19.1, mice 3.19.0, ggplot2 4.0.3.
All random operations are seeded, and results do not depend on the number of cores used by `04`. Tree
building with a different xgboost version can change the last decimals of the importance values.

## Folders

* `output/` results of the last complete run (CSV; two jpeg figures)
* `expected/` the Supplementary Tables as printed, used by script 06
* `cache/` intermediate objects (can be deleted; they are rebuilt)
