# -*- coding: utf-8 -*-
"""
Step 03: Integrated Cross-Figure Phylodynamic & Ecological Statistics.

Merges 14-sublineage spatial breadth (Step 01), post-QC BEAST/BSP metrics (Step 02),
and sublineage drug-resistance burdens to produce:
  - Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv
  - Table_S2b_CrossFigure_Statistical_Tests.tsv
    (including all Spearman/Pearson correlations, Global vs Non-Global Mann-Whitney U
     group contrasts, and complete 14-sublineage Jackknife leave-one-out sensitivity tests)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, ensure_output_dirs


def run_cross_figure_statistics():
    print("[01_analysis | Step 3] Integrating spatial, BSP, BEAST & DR tables and running cross-figure tests...")
    ensure_output_dirs()

    df_geo = pd.read_csv(TABLES_DIR / "Table_S2_spatial_metrics_L2.tsv", sep="\t")
    df_bsp = pd.read_csv(TABLES_DIR / "Table_S2_bsp_metrics.tsv", sep="\t")
    df_tree = pd.read_csv(TABLES_DIR / "Table_S4_BEAST_Tree_QC_and_Outlier_Audit.tsv", sep="\t")
    df_dr = pd.read_csv(
        DATA_DIR / "Extended Data Fig. 4" / "sublineage_resistance_burden.tsv", sep="\t"
    ).rename(columns={"n_samples": "DR_n_samples"})

    df_m = (
        df_geo.merge(df_bsp, on="sublineage")
        .merge(
            df_tree[["sublineage", "BEAST_ntax", "root_height_yr", "root_top1_loc", "root_top1_PP", "QC_status"]],
            on="sublineage",
        )
        .merge(df_dr, on="sublineage")
    )
    df_m.to_csv(TABLES_DIR / "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv", sep="\t", index=False)

    stat_rows = []
    # 1. Pairwise Spearman & Pearson correlations (full cohort n = 14 + leave-L4.1-out subset)
    for label, sub_df in [
        ("All_14_Sublineages", df_m),
        ("Sensitivity_Excl_L4.1", df_m[df_m["sublineage"] != "lineage4.1"]),
    ]:
        for xcol in ["n_countries", "shannon_H", "top3_share", "any_DR_rate", "MDR_rate", "pre_XDR_rate"]:
            for ycol in ["retain_vs_peak", "retain_vs_1950", "retain_vs_2000"]:
                rho, p_sp = stats.spearmanr(sub_df[xcol], sub_df[ycol])
                r, p_pe = stats.pearsonr(sub_df[xcol], sub_df[ycol])
                stat_rows.append({
                    "Analysis_Set": label,
                    "Test_Type": "Correlation",
                    "Variable_X": xcol,
                    "Variable_Y": ycol,
                    "Statistic_Primary": f"Spearman_rho={rho:.4f}",
                    "P_Value_Primary": round(float(p_sp), 6),
                    "Statistic_Secondary": f"Pearson_r={r:.4f}",
                    "P_Value_Secondary": round(float(p_pe), 6),
                    "Significant_0.05": "YES" if (p_sp < 0.05 or p_pe < 0.05) else "NO",
                })

    # 2. Two-sided Mann-Whitney U contrasts (Global n = 5 vs Non-Global n = 9)
    g_df = df_m[df_m["ecological_group"] == "Global"]
    ng_df = df_m[df_m["ecological_group"] != "Global"]
    r_df = df_m[df_m["ecological_group"] == "Rooted"]
    for col in [
        "n_countries", "shannon_H", "top3_share",
        "retain_vs_peak", "retain_vs_1950", "retain_vs_2000",
        "any_DR_rate", "MDR_rate", "pre_XDR_rate",
    ]:
        u, p = stats.mannwhitneyu(g_df[col], ng_df[col], alternative="two-sided")
        stat_rows.append({
            "Analysis_Set": "Global(n=5)_vs_NonGlobal(n=9)",
            "Test_Type": "MannWhitneyU_TwoSided",
            "Variable_X": "ecological_group",
            "Variable_Y": col,
            "Statistic_Primary": (
                f"Global_med={g_df[col].median():.4f}; "
                f"NonGlobal_med={ng_df[col].median():.4f}; "
                f"Rooted_med={r_df[col].median():.4f}"
            ),
            "P_Value_Primary": round(float(p), 6),
            "Statistic_Secondary": f"U={u:.1f}",
            "P_Value_Secondary": np.nan,
            "Significant_0.05": "YES" if p < 0.05 else "NO",
        })

    # 3. Complete 14-sublineage Jackknife Leave-One-Out (n = 13 per iteration)
    for s_out in df_m["sublineage"]:
        sub = df_m[df_m["sublineage"] != s_out]
        rho, p_sp = stats.spearmanr(sub["n_countries"], sub["retain_vs_2000"])
        r, p_pe = stats.pearsonr(sub["n_countries"], sub["retain_vs_2000"])
        g_i = sub[sub["ecological_group"] == "Global"]["retain_vs_2000"]
        ng_i = sub[sub["ecological_group"] != "Global"]["retain_vs_2000"]
        _, p_mwu = stats.mannwhitneyu(g_i, ng_i, alternative="two-sided")
        lbl = s_out.replace("lineage", "L")
        stat_rows.append({
            "Analysis_Set": f"LeaveOneOut_Excl_{lbl}",
            "Test_Type": "Jackknife_LOO_n13",
            "Variable_X": "n_countries",
            "Variable_Y": "retain_vs_2000",
            "Statistic_Primary": f"Spearman_rho={rho:.4f}; Pearson_r={r:.4f}",
            "P_Value_Primary": round(float(p_sp), 6),
            "Statistic_Secondary": f"MWU_Global_vs_NonGlobal_P={p_mwu:.6f}",
            "P_Value_Secondary": round(float(p_mwu), 6),
            "Significant_0.05": "YES" if (p_sp < 0.05 and p_mwu < 0.05) else "NO",
        })

    df_stats = pd.DataFrame(stat_rows)
    df_stats.to_csv(TABLES_DIR / "Table_S2b_CrossFigure_Statistical_Tests.tsv", sep="\t", index=False)

    # Print key headline statistics for quick verification
    rho_main, p_main = stats.spearmanr(df_m["n_countries"], df_m["retain_vs_2000"])
    _, p_mwu_main = stats.mannwhitneyu(g_df["retain_vs_2000"], ng_df["retain_vs_2000"], alternative="two-sided")
    print(f"  Primary test: breadth vs 21st-c. retention Spearman rho={rho_main:+.3f} (P={p_main:.4f}); MWU P={p_mwu_main:.4f}")
    print(f"  Saved Table_S2 ({len(df_m)} sublineages) & Table_S2b ({len(df_stats)} tests, incl. 14 LOO jackknife sets).")
    return df_m, df_stats


def main():
    run_cross_figure_statistics()


if __name__ == "__main__":
    main()
