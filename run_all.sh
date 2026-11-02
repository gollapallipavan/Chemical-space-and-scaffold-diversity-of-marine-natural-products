#!/usr/bin/env bash
# Reproduces Steps 3-5 from the shipped supplementary data (no SDF needed).
# Steps 1-2 need the raw CMNPD SDF; see README.
set -euo pipefail
cd "$(dirname "$0")/scripts"
python 03_scaffold_mcps.py   | tee ../results/log_03_scaffold_mcps.txt
python 04_extended_analysis.py 2>&1 | grep -v -i deprecat | tee ../results/log_04_extended_analysis.txt
python 05_make_figures.py    | tee ../results/log_05_make_figures.txt
