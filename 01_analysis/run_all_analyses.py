# -*- coding: utf-8 -*-
"""
Master Analysis Pipeline Runner (`01_analysis/run_all_analyses.py`)

Executes Steps 01 through 05 in sequence:
  01. Cohort filtering funnel (162,869 -> 6,712) & L2/L3 spatial niche metrics
  02. BEAST MCC tree QC, post-QC L4.1 outlier pruning (2405024009) & BSP alignment
  03. Cross-figure spatial-vs-phylodynamic correlations, MWU & 14-fold Jackknife LOO
  04. Convergent selection scan & 18-gene clonal clade-hitchhiking filtration
  05. Sublineage drug-resistance burden & low-fitness-cost mutation enrichment

Usage:
  export MTBC_DATA_DIR=/path/to/raw_data
  python 01_analysis/run_all_analyses.py
"""

import importlib
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
sys.path.insert(0, str(THIS_DIR.parent))

from config import DATA_DIR, RESULTS_DIR


STEPS = [
    "01_cohort_audit_and_spatial_metrics",
    "02_beast_qc_and_L41_outlier_pruning",
    "03_cross_figure_phylodynamic_statistics",
    "04_convergent_selection_and_hitchhiking_filter",
    "05_drug_resistance_enrichment_analysis",
]


def main():
    print("=" * 72)
    print("MTBC Global Evolution — Complete Analysis Pipeline (01_analysis)")
    print(f"  Input data directory : {DATA_DIR}")
    print(f"  Output results root  : {RESULTS_DIR}")
    print("=" * 72)

    for mod_name in STEPS:
        mod = importlib.import_module(mod_name)
        mod.main()

    print("=" * 72)
    print("ALL ANALYSIS STEPS COMPLETED SUCCESSFULLY.")
    print("=" * 72)


if __name__ == "__main__":
    main()
