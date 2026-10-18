# PAnACEA-H: analysis code

Code accompanying:

> Development and internal validation of the PAnACEA-H score for predicting
> hospital admission at Emergency Department triage in adult patients
> presenting to United States hospitals: a machine learning analysis
> (*Frontiers in Public Health*, under review)

## Repository contents

| Path | Purpose |
|---|---|
| `panacea_h_reproducibility.R` | Full pipeline: feature ranking, score derivation, 10-fold CV, held-out validation, bootstrap CIs, decision curve analysis |
| `regenerate_figures.R` | Main-text Figures 1-3 and Supplementary Figure 3 (see "Scope" below) |
| `data_for_figures/` | Small precomputed inputs for the figure script |
| `reproducibility_output/` | Reference output of a complete run of the pipeline (for comparison) |

## Scope: what this repository does and does not contain

**Included**

- Feature ranking by information gain (Figure 2; Supplementary Table 3)
- Pooled 5-fold cross-validated performance of the XGBoost ranking model
  (Supplementary Table 12)
- Score item weights and cutoffs (Table 2)
- Internal validation by 10-fold cross-validation on the training set
  (Supplementary Table 15)
- Held-out test-set validation: threshold-level performance (Table 3), the
  ROC curve (Figure 3), and bootstrap confidence intervals (Supplementary
  Table 16, Supplementary Figure 2)
- Decision curve analysis of the score against Admit All / Admit None
  (the pipeline's own version; Supplementary Figure 4 in the manuscript also
  shows ESI alone)
- Figures 1-3 and Supplementary Figure 3 (calibration) as they appear in the
  revised manuscript

The output files written by the pipeline keep the supplementary numbering of
the original submission in their names (e.g. `supplementary_table4_...`); the
numbers above are those of the revised manuscript.

**Not included**

- Table 1 (descriptive baseline characteristics): purely descriptive.
- The additional analyses added during peer review (calibration computation,
  subgroup and site-level performance, optimism correction, selection
  stability, comparator models, and the remaining new supplementary tables).
  Only the calibration data needed to draw Supplementary Figure 3 are provided
  (`data_for_figures/R2_calibration_*.csv`), as precomputed files; the code
  that generated them is not part of this release.
- Figure 1 is drawn from standardized mean differences transcribed from
  Supplementary Tables 9-11; they are not recomputed from the raw data.

The ROC curve in main-text Figure 3 is drawn by `regenerate_figures.R` from the
threshold table written by the pipeline (the pipeline's own `figure3_roc_test_set.jpeg`
is an equivalent, unformatted version).

## Data

The dataset is not redistributed here. Download it from the source study's
repository:

- Hong WS, Haimovich AD, Taylor RA (2018). Predicting hospital admission at
  emergency department triage using machine learning. *PLoS ONE* 13(7):
  e0201016. <https://doi.org/10.1371/journal.pone.0201016>
- Dataset: <https://github.com/yaleemmlc/admissionprediction> ->
  `Results/5v_cleandf.RData` (560,486 visits x 972 variables)

Place `5v_cleandf.RData` in the repository root (or edit `DATA_FILE` at the
top of the script).

## How to run

From the repository root:

```r
install.packages(c("dplyr", "tidyr", "caret", "xgboost", "pROC",
                   "mice", "ggplot2", "dcurves"))
source("panacea_h_reproducibility.R")   # tens of minutes (multiple imputation)
source("regenerate_figures.R")          # seconds; writes ./figures/
```

The pipeline writes all output to `./reproducibility_output/`, overwriting
the reference files shipped here. The figure script does not need the raw
dataset.

## A note on reproducibility

The `max_depth` grid search included for transparency in Section 4b of the
pipeline selects `max_depth = 8` (cross-validated AUC ~0.92) as its own best
choice, not the `max_depth = 3` stated in the manuscript's Methods and fixed
explicitly in this script for the published ranking model. Using
`max_depth = 3` reproduces the manuscript's reported pooled AUC of 91.4%.
The diagnostic grid search is in
`reproducibility_output/diagnostic_max_depth_grid_search.csv`.

## Software

Developed and tested under R 4.2.3. Package versions are not pinned. The
`mice` step uses random draws (predictive mean matching), so re-running may
reproduce the reported figures to within a few hundredths of a percentage
point rather than bit-for-bit, depending on package versions.

## License

Code: MIT (see `LICENSE`). The underlying dataset is subject to the terms of
its original source (Hong et al., 2018).

## Citation

See `CITATION.cff`. Version 1.0.0 is archived on Zenodo:
[10.5281/zenodo.23450955](https://doi.org/10.5281/zenodo.23450955)
(concept DOI covering all versions: [10.5281/zenodo.23450954](https://doi.org/10.5281/zenodo.23450954)).
