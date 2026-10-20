#!/usr/bin/env bash
# Runs the analyses added during peer review, in order, from the repository root:
#   PANACEA_DATA=/path/to/5v_cleandf.RData bash revision_analyses/run_all.sh
# The whole folder takes about 20 minutes on an 18-core machine (04 runs in parallel).
set -euo pipefail
for s in 01_mapping_calibration_bands 02_comparators 03_subgroups 04_optimism_and_selection 05_sensitivity_analyses 06_verify_against_manuscript; do
  echo "=== $s ($(date '+%H:%M')) ==="
  Rscript "revision_analyses/$s.R"
done
