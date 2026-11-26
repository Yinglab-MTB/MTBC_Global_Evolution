# -*- coding: utf-8 -*-
"""
Step 05: Sublineage Drug-Resistance Burden & Low-Fitness-Cost Mutation Enrichment.

Summarizes drug-resistance phenotypes across the 14-sublineage cohort (n = 23,533),
verifies two-sided Fisher's exact tests and 95% Woolf confidence intervals for the
canonical low-fitness-cost resistance mutations enriched in Global (n = 15,257) vs.
Rooted (n = 941) clades, and exports clean summary tables:
  - Table_S5a_Sublineage_DR_Burden.tsv
  - Table_S5b_Global_vs_Rooted_LowCost_DR_Enrichment.tsv
  - Table_S5c_Within_Lineage_DR_Epistasis_Summary.tsv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, GLOBAL_SET, ROOTED_SET, ensure_output_dirs


def run_dr_enrichment_analysis():
    print("[01_analysis | Step 5] Summarizing sublineage DR burdens & low-fitness-cost mutation enrichments...")
    ensure_output_dirs()

    dr_dir = DATA_DIR / "Extended Data Fig. 4"
    df_bur = pd.read_csv(dr_dir / "sublineage_resistance_burden.tsv", sep="\t")
    df_drg = pd.read_csv(dr_dir / "DR_mutation_global_enriched.tsv", sep="\t")
    df_wl = pd.read_csv(dr_dir / "DR_mutation_within_lineage_target_enriched.tsv", sep="\t")

    # 1. Sublineage DR burden table with ecological group labels
    df_bur["ecological_group"] = df_bur["sublineage"].map(
        lambda s: "Global" if s in GLOBAL_SET else ("Rooted" if s in ROOTED_SET else "Intermediate")
    )
    df_bur = df_bur.sort_values("MDR_rate", ascending=False).reset_index(drop=True)
    df_bur.to_csv(TABLES_DIR / "Table_S5a_Sublineage_DR_Burden.tsv", sep="\t", index=False)

    # 2. Verify Global (n = 15,257) vs. Rooted (n = 941) Fisher exact tests & 95% OR CIs
    n_global = int(df_bur[df_bur["ecological_group"] == "Global"]["n_samples"].sum())
    n_rooted = int(df_bur[df_bur["ecological_group"] == "Rooted"]["n_samples"].sum())

    ci_lo, ci_hi, p_ver = [], [], []
    for _, r in df_drg.iterrows():
        a = int(r["global"])
        b = n_global - a
        c = int(r["rooted"])
        d = n_rooted - c
        or_val, p_val = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
        se_log_or = np.sqrt(1.0 / max(a, 0.5) + 1.0 / max(b, 0.5) + 1.0 / max(c, 0.5) + 1.0 / max(d, 0.5))
        ci_lo.append(round(float(np.exp(np.log(or_val) - 1.96 * se_log_or)), 4))
        ci_hi.append(round(float(np.exp(np.log(or_val) + 1.96 * se_log_or)), 4))
        p_ver.append(float(p_val))

    df_drg_out = df_drg.copy()
    df_drg_out["OR_95CI_low"] = ci_lo
    df_drg_out["OR_95CI_high"] = ci_hi
    df_drg_out["fisher_p_verified"] = p_ver
    df_drg_out = df_drg_out.sort_values("FDR", ascending=True).reset_index(drop=True)
    df_drg_out.to_csv(TABLES_DIR / "Table_S5b_Global_vs_Rooted_LowCost_DR_Enrichment.tsv", sep="\t", index=False)

    # 3. Within-lineage DR enrichment summary (L2.2 vs L2.1 & global L4 sublineages)
    df_wl.to_csv(TABLES_DIR / "Table_S5c_Within_Lineage_DR_Epistasis_Summary.tsv", sep="\t", index=False)

    # Report Global vs Non-Global MDR-TB Mann-Whitney U test
    g_mdr = df_bur[df_bur["ecological_group"] == "Global"]["MDR_rate"]
    ng_mdr = df_bur[df_bur["ecological_group"] != "Global"]["MDR_rate"]
    _, p_mdr = stats.mannwhitneyu(g_mdr, ng_mdr, alternative="two-sided")
    comp_counts = df_wl["comparison"].value_counts().to_dict() if "comparison" in df_wl.columns else {}
    print(f"  Global vs Non-Global MDR-TB median: {g_mdr.median()*100:.1f}% vs {ng_mdr.median()*100:.1f}% (MWU P={p_mdr:.3f}, n.s.)")
    print(f"  Saved Table_S5a ({len(df_bur)} clades), Table_S5b ({len(df_drg_out)} low-cost mutations), Table_S5c ({len(df_wl)} within-lineage variants; {comp_counts}).")
    return df_bur, df_drg_out, df_wl


def main():
    run_dr_enrichment_analysis()


if __name__ == "__main__":
    main()
