# -*- coding: utf-8 -*-
"""
Master Figure Generation Runner (`02_plotting/run_all_figures.py`)

Generates all 10 publication-grade figures (vector PDF + 400 DPI PNG):
  - Fig. 1 : Global MTBC Atlas (Robinson donut map, dual-ring sunburst, trajectory, Marimekko)
  - Fig. 2 : Second-Level Sublineage Geography (14 Kriging contour maps + breadth summary)
  - Fig. 3 : Third-Level Nested Specialization (6 L3 Kriging maps, 30-clade polar ring, dumbbell)
  - Fig. 4 : Demographic Resilience & BSP Dynamics (overlay, 14 sparklines, regression, rainclouds)
  - Fig. 5 : Sublineage Phylogeography Trees (6 circular BEAST trees, tMRCA caterpillar, root PP)
  - Fig. 6 : Convergent Adaptation & Low-Cost DR (4.41 Mb Circos, Fst matrix, OR forest, 2-phase model)
  - Extended Data Figs. 1–4 : Cohort funnel, ecological niche metrics, post-QC/jackknife/artifacts, DR architecture

Prerequisite:
  Run `python 01_analysis/run_all_analyses.py` first so summary tables and
  post-QC L4.1 BEAST files exist under `./results/`.

Usage:
  export MTBC_DATA_DIR=/path/to/raw_data
  python 02_plotting/run_all_figures.py
"""

import importlib
import sys
from pathlib import Path

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
sys.path.insert(0, str(THIS_DIR.parent))

from config import DATA_DIR, FIGURES_DIR


FIGURE_SCRIPTS = [
    "plot_fig1_global_atlas",
    "plot_fig2_second_level_geography",
    "plot_fig3_third_level_nested_specialization",
    "plot_fig4_demographic_resilience_bsp",
    "plot_fig5_phylogeography_trees",
    "plot_fig6_convergent_adaptation_and_dr",
    "plot_extended_data_figs",
]


def main():
    print("=" * 72)
    print("MTBC Global Evolution — Complete Figure Generation Suite (02_plotting)")
    print(f"  Input data directory : {DATA_DIR}")
    print(f"  Output figures dir   : {FIGURES_DIR}")
    print("=" * 72)

    for mod_name in FIGURE_SCRIPTS:
        mod = importlib.import_module(mod_name)
        mod.main()

    print("=" * 72)
    print(f"ALL 10 PUBLICATION FIGURES SAVED TO: {FIGURES_DIR}")
    print("=" * 72)


if __name__ == "__main__":
    main()
