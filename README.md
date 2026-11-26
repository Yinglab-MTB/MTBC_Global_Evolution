# Global Sublineage Evolution and Antibiotic-Era Demographic Resilience of *Mycobacterium tuberculosis*

Source code accompanying the manuscript:
> **Convergent cell-envelope adaptation and low-cost drug resistance underpin the global success and modern demographic resilience of *Mycobacterium tuberculosis* sublineages**

This repository is organized into two clearly separated modules for transparent peer review and full reproducibility:
1. **`01_analysis/` — Statistical, Phylodynamic QC & Population-Genomic Analysis Code（分析代码）**
2. **`02_plotting/` — Publication-Grade Figure Generation Code（绘图代码）**

---

## Repository Structure

```text
GitHub_MTBC_Global_Evolution/
├── README.md                                          # Documentation & reproduction guide
├── LICENSE                                            # MIT License
├── requirements.txt                                   # Python dependencies
├── config.py                                          # Portable path configuration (env-overridable)
├── assets/
│   └── world_map_polygons.csv                         # Coastline geometry for Robinson projection maps
│
├── 01_analysis/                                       # [Category 1] Data analysis & statistical pipelines
│   ├── 01_cohort_audit_and_spatial_metrics.py         # 4-stage cohort filtering (162,869 → 6,712) & L2/L3 spatial metrics
│   ├── 01b_prepare_surfer_kriging_inputs.R            # R script exporting lat/lon/frequency tables for Surfer Kriging
│   ├── 02_beast_qc_and_L41_outlier_pruning.py         # BEAST MCC tree QC (H1/H2, tMRCA, ancestral PP), L4.1 outlier pruning (2405024009) & FASTA filter CLI
│   ├── 03_cross_figure_phylodynamic_statistics.py     # Breadth vs. Ne retention correlations, MWU tests & 14-sublineage Jackknife LOO sensitivity
│   ├── 04_convergent_selection_and_hitchhiking_filter.py  # Multi-evidence selection scan (Δf, Fst, π, Tajima's D) & 18-gene clonal hitchhiking filter
│   ├── 05_drug_resistance_enrichment_analysis.py      # 14-sublineage DR burden, low-fitness-cost mutation enrichment & within-L2/L4 epistasis
│   └── run_all_analyses.py                            # Master script executing all analysis steps sequentially
│
└── 02_plotting/                                       # [Category 2] Publication figure generation (PDF + 400 DPI PNG)
    ├── nature_v2_style.py                             # Unified Nature editorial palette, Robinson projection, Kriging & circular BEAST tree engine
    ├── plot_fig1_global_atlas.py                      # Main Fig. 1 (Global donut map, dual-ring sunburst, temporal trajectory, Marimekko mosaic)
    ├── plot_fig2_second_level_geography.py            # Main Fig. 2 (14 second-level sublineage Kriging contour maps & breadth summary)
    ├── plot_fig3_third_level_nested_specialization.py # Main Fig. 3 (6 contrasting L3 Kriging maps, 30-clade polar ring & nested specialization)
    ├── plot_fig4_demographic_resilience_bsp.py        # Main Fig. 4 (Normalized BSP overlay, 14 BSP sparklines, regression & rainclouds)
    ├── plot_fig5_phylogeography_trees.py              # Main Fig. 5 (6 circular BEAST MCC trees, tMRCA 95% HPD caterpillar & ancestral origin bars)
    ├── plot_fig6_convergent_adaptation_and_dr.py      # Main Fig. 6 (4.41 Mb Circos scan, Tier 1/2 Fst matrix, low-cost DR forest & two-phase model)
    ├── plot_extended_data_figs.py                     # Extended Data Figs. 1–4 (Funnel, niche quadrant, QC & 14-set jackknife, DR architecture)
    └── run_all_figures.py                             # Master script rendering all 10 publication figures
```

---

## Mapping of Code to Manuscript Outputs

| Manuscript Output | Analysis Script (`01_analysis/`) | Plotting Script (`02_plotting/`) |
| :--- | :--- | :--- |
| **Fig. 1 & Extended Data Fig. 1** (`Supplementary Table 1/1b`) | `01_cohort_audit_and_spatial_metrics.py` | `plot_fig1_global_atlas.py`, `plot_extended_data_figs.py` |
| **Fig. 2, Fig. 3 & Extended Data Fig. 2** (`Supplementary Table 2/2c`) | `01_cohort_audit_and_spatial_metrics.py`, `01b_prepare_surfer_kriging_inputs.R` | `plot_fig2_second_level_geography.py`, `plot_fig3_third_level_nested_specialization.py`, `plot_extended_data_figs.py` |
| **Fig. 4 & Extended Data Fig. 3b** (`Supplementary Table 2/2b`) | `02_beast_qc_and_L41_outlier_pruning.py`, `03_cross_figure_phylodynamic_statistics.py` | `plot_fig4_demographic_resilience_bsp.py`, `plot_extended_data_figs.py` |
| **Fig. 5 & Extended Data Fig. 3a** (`Supplementary Table 4`) | `02_beast_qc_and_L41_outlier_pruning.py` | `plot_fig5_phylogeography_trees.py`, `plot_extended_data_figs.py` |
| **Fig. 6 & Extended Data Fig. 3c** (`Supplementary Table 3/3b`) | `04_convergent_selection_and_hitchhiking_filter.py`, `05_drug_resistance_enrichment_analysis.py` | `plot_fig6_convergent_adaptation_and_dr.py`, `plot_extended_data_figs.py` |
| **Extended Data Fig. 4** (`Supplementary Table 5`) | `05_drug_resistance_enrichment_analysis.py` | `plot_extended_data_figs.py` |

---

## Installation

### Python Environment (Tested on Python 3.9–3.11)
```bash
git clone https://github.com/<your-username>/MTBC-Global-Evolution.git
cd MTBC-Global-Evolution
pip install -r requirements.txt
```

### Optional R Dependencies (for `01_analysis/01b_prepare_surfer_kriging_inputs.R`)
```r
install.packages(c("readxl", "openxlsx", "dplyr"))
```

---

## Input Data Layout

Place the study source data folders under `./data/` (or point `MTBC_DATA_DIR` to your existing data directory):

```text
data/
├── Fig. 1/                  # data.xlsx, data.csv
├── Fig. 2/                  # 二级谱系分布频率.xlsx (country-level L2/L3 frequency sheets)
├── Fig. 3/                  # 二级谱系分布频率.xlsx (3rd-level sublineage frequency sheets)
├── Fig. 4/                  # *.BSP.txt (Tracer Bayesian Skyline tables for 14 sublineages)
├── Fig. 5/                  # *_location_tree_with_trait.tree (BEAST MCC discrete-trait trees)
├── Fig. 6/                  # candidate_gene_*.tsv (Fst, pi, Tajima's D, coordinates, summaries)
└── Extended Data Fig. 4/    # sublineage_resistance_burden.tsv, DR_mutation_*_enriched.tsv
```

---

## Quick Start（快速运行）

### 1. Run All Analyses（运行全部分析代码）
```bash
# Linux / macOS
export MTBC_DATA_DIR=/path/to/raw_data
python 01_analysis/run_all_analyses.py

# Windows PowerShell
$env:MTBC_DATA_DIR = "C:\path\to\raw_data"
python 01_analysis/run_all_analyses.py
```
Outputs are saved to `./results/tables/` (`Table_S1` through `Table_S5` TSVs) and `./results/cleaned_L41_beast/` (post-QC 499-tip `L4.1` BEAST tree & crown-aligned BSP after removing basal outlier `2405024009` exclusively from BEAST analyses).

### 2. Generate All 10 Publication Figures（运行全部绘图代码）
```bash
python 02_plotting/run_all_figures.py
```
All 10 figures (Main **Figs. 1–6** and **Extended Data Figs. 1–4**) are exported as both vector `.pdf` and 400 DPI `.png` into `./results/figures/`.

### 3. Optional: Pre-BEAST FASTA / Trait Filtering for `L4.1` Outlier `2405024009`
```bash
python 01_analysis/02_beast_qc_and_L41_outlier_pruning.py \
    --filter-fasta L4.1_aln.fasta --out-fasta L4.1_aln_clean499.fasta \
    --filter-meta L4.1_traits.tsv --out-meta L4.1_traits_clean499.tsv
```

---

## Key Methodological Notes
- **Strict Post-QC Scope**: Isolate `2405024009` (anomalous basal long branch, pre-QC $H_1/H_2 = 2.06$) is excluded **exclusively** from `L4.1` BEAST phylodynamic and phylogeographic inference (`L4.1` $n = 499$, $H_1/H_2 = 1.02$; Stage 4 total $n = 6,712$), while all non-BEAST cohorts (`162,869` global atlas, `107,309` geolocated atlas, `23,533` comparative selection & resistance cohort) remain unchanged.
- **Clonal Hitchhiking Control**: `04_convergent_selection_and_hitchhiking_filter.py` explicitly identifies and filters 18 single-SNP ancestral clade-hitchhiking artifacts (Cluster A `L4.7–L4.9/H37Rv`, Cluster B `L4.6–L4.9`, Cluster C `L4`-dominant) prior to calling the 9 Tier-1 cross-lineage (`L2 + L4`) multi-hit convergent adaptive genes.

## License
Released under the [MIT License](LICENSE).
