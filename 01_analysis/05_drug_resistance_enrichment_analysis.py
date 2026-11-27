# -*- coding: utf-8 -*-
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, GLOBAL_SET, ROOTED_SET, ensure_output_dirs


def run_dr_enrichment_analysis():
    ensure_output_dirs()

    dr_dir = DATA_DIR
    df_bur = pd.read_csv(dr_dir / "sublineage_drug_resistance_burden.tsv", sep="\t")
    df_drg = pd.read_csv(dr_dir / "drug_resistance_mutations_global_vs_rooted_enriched.tsv", sep="\t")
    df_wl = pd.read_csv(dr_dir / "drug_resistance_mutations_within_lineage_enriched.tsv", sep="\t")

    df_bur["ecological_group"] = df_bur["sublineage"].map(
        lambda s: "Global" if s in GLOBAL_SET else ("Rooted" if s in ROOTED_SET else "Intermediate")
    )
    df_bur = df_bur.sort_values("MDR_rate", ascending=False).reset_index(drop=True)
    df_bur.to_csv(TABLES_DIR / "Table_S5a_Sublineage_DR_Burden.tsv", sep="\t", index=False)

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

    df_wl.to_csv(TABLES_DIR / "Table_S5c_Within_Lineage_DR_Epistasis_Summary.tsv", sep="\t", index=False)

    g_mdr = df_bur[df_bur["ecological_group"] == "Global"]["MDR_rate"]
    ng_mdr = df_bur[df_bur["ecological_group"] != "Global"]["MDR_rate"]
    _, p_mdr = stats.mannwhitneyu(g_mdr, ng_mdr, alternative="two-sided")
    comp_counts = df_wl["comparison"].value_counts().to_dict() if "comparison" in df_wl.columns else {}
    return df_bur, df_drg_out, df_wl


def main():
    run_dr_enrichment_analysis()


if __name__ == "__main__":
    main()
