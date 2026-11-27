# **Convergent cell-envelope adaptation and low-cost drug resistance underpin the global success and modern demographic resilience of *Mycobacterium tuberculosis* sublineages**

Analysis code accompanying the manuscript *Convergent cell-envelope adaptation and low-cost drug resistance underpin the global success and modern demographic resilience of Mycobacterium tuberculosis sublineages*.

## Getting started

From the repository root, install the Python dependencies and run the analyses after placing the study inputs in `data/`.

```sh
python -m pip install -r requirements.txt
python 01_analysis/run_all_analyses.py
```

Set `MTBC_DATA_DIR` and `MTBC_RESULTS_DIR` to use a different input directory or output directory. Relative paths supplied through these variables are resolved from the current working directory.

```sh
# Linux or macOS
export MTBC_DATA_DIR="/path/to/study-data"
export MTBC_RESULTS_DIR="/path/to/results"
python 01_analysis/run_all_analyses.py
```

```powershell
# Windows PowerShell
$env:MTBC_DATA_DIR = 'C:\path\to\study-data'
$env:MTBC_RESULTS_DIR = 'C:\path\to\results'
python 01_analysis/run_all_analyses.py
```

## Study inputs

The `data/` directory is excluded from Git. Supply these files in `data/` or in the directory configured by `MTBC_DATA_DIR`.

```text
data/
├── MTBC_global_cohort_metadata_and_summaries.xlsx
├── MTBC_country_lineage_counts_coordinates.csv
├── MTBC_sublineage_country_frequencies_levels_2_3.xlsx
├── sublineage_drug_resistance_burden.tsv
├── selection_candidate_genes_summary_all_comparisons.tsv
├── selection_candidate_genes_H37Rv_coordinates.tsv
├── selection_candidate_genes_Fst_pi_TajimaD_all_contrasts.tsv
├── drug_resistance_mutations_global_vs_rooted_enriched.tsv
├── drug_resistance_mutations_within_lineage_enriched.tsv
├── phylogeographic_trees/
│   └── lineage*_location_tree_with_trait.tree
└── bayesian_skyline/
    └── L*.BSP.txt
```

## Optional Surfer input export

Install the R packages in an R session.

```r
install.packages(c("readxl", "openxlsx", "dplyr"))
```

Run the export separately from the repository root.

```sh
Rscript 01_analysis/01b_prepare_surfer_kriging_inputs.R
```

## License

The code is released under the [MIT License](LICENSE).
