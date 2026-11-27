# -*- coding: utf-8 -*-
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, SUB14, ensure_output_dirs


# Full-cohort correlations, group contrasts, and leave-one-out tests.
X_COLS = ["n_countries", "shannon_H", "top3_share", "any_DR_rate", "MDR_rate", "pre_XDR_rate"]
Y_COLS = ["retain_vs_peak", "retain_vs_1950", "retain_vs_2000"]


def validate_inputs(frames):
    expected = set(SUB14)
    if len(expected) != 14 or "lineage4.1" not in expected:
        raise ValueError("SUB14 must contain 14 distinct sublineages including lineage4.1")
    for name, frame in frames.items():
        ids = frame["sublineage"]
        if ids.isna().any() or ids.duplicated().any():
            raise ValueError(f"{name}: missing or duplicate sublineage IDs")
        missing = expected - set(ids)
        if missing:
            raise ValueError(f"{name}: missing sublineages: {sorted(missing)}")
    tree = frames["tree"]
    n = tree.loc[tree["sublineage"] == "lineage4.1", "BEAST_ntax"].iloc[0]
    if n != 499:
        raise ValueError(f"L4.1 expected 499 tree taxa, found {n}")
    # ntax alone cannot establish that the correct sample was removed,
    # nor that BSP was inferred from the same 499 taxa.


def correlations(frame, xcol, ycol):
    x, y = frame[xcol], frame[ycol]
    if len(frame) < 3 or x.nunique() < 2 or y.nunique() < 2:
        raise ValueError(f"Undefined correlation: {xcol} vs {ycol}; n={len(frame)}")
    rho, p_sp = stats.spearmanr(x, y)
    r, p_pe = stats.pearsonr(x, y)
    return rho, p_sp, r, p_pe


def bh_adjust(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(p)
    out[order] = np.minimum(1.0, np.minimum.accumulate(ranked[::-1])[::-1])
    return out


def run_cross_figure_statistics():
    ensure_output_dirs()

    df_geo = pd.read_csv(TABLES_DIR / "Table_S2_spatial_metrics_L2.tsv", sep="\t")
    df_bsp = pd.read_csv(TABLES_DIR / "Table_S2_bsp_metrics.tsv", sep="\t")
    df_tree = pd.read_csv(TABLES_DIR / "Table_S4_BEAST_Tree_QC_and_Outlier_Audit.tsv", sep="\t")
    df_dr = pd.read_csv(
        DATA_DIR / "sublineage_drug_resistance_burden.tsv", sep="\t"
    ).rename(columns={"n_samples": "DR_n_samples"})

    frames = {"geo": df_geo, "bsp": df_bsp, "tree": df_tree, "dr": df_dr}
    validate_inputs(frames)
    df_geo, df_bsp, df_tree, df_dr = (
        f.loc[f["sublineage"].isin(SUB14)].copy() for f in frames.values()
    )

    df_m = (
        df_geo.merge(df_bsp, on="sublineage", validate="one_to_one")
        .merge(
            df_tree[["sublineage", "BEAST_ntax", "root_height_yr", "root_top1_loc", "root_top1_PP", "QC_status"]],
            on="sublineage", validate="one_to_one",
        )
        .merge(df_dr, on="sublineage", validate="one_to_one")
    )
    if len(df_m) != 14:
        raise ValueError(f"Expected 14 merged sublineages, found {len(df_m)}")
    for col in X_COLS + Y_COLS:
        df_m[col] = pd.to_numeric(df_m[col], errors="raise")
        if not np.isfinite(df_m[col].to_numpy(dtype=float)).all():
            raise ValueError(f"Missing or non-finite values in {col}")
    if df_m["ecological_group"].isna().any():
        raise ValueError("Missing ecological_group")
    if (df_m["ecological_group"] == "Global").sum() != 5 or (df_m["ecological_group"] == "Rooted").sum() == 0:
        raise ValueError("Unexpected ecological group composition; verify upstream classifications")

    stat_rows = []
    # 1. Pairwise Spearman & Pearson correlations (full cohort n = 14)
    for label, sub_df in [
        ("All_14_Sublineages", df_m),
    ]:
        for xcol in X_COLS:
            for ycol in Y_COLS:
                rho, p_sp, r, p_pe = correlations(sub_df, xcol, ycol)
                stat_rows.append({
                    "Analysis_Set": label,
                    "Test_Type": "Correlation",
                    "Variable_X": xcol,
                    "Variable_Y": ycol,
                    "Statistic_Primary": f"Spearman_rho={rho:.4f}",
                    "P_Value_Primary": float(p_sp),
                    "Statistic_Secondary": f"Pearson_r={r:.4f}",
                    "P_Value_Secondary": float(p_pe),
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
            "P_Value_Primary": float(p),
            "Statistic_Secondary": f"U={u:.1f}",
            "P_Value_Secondary": np.nan,
            "Significant_0.05": "YES" if p < 0.05 else "NO",
        })

    # 3. Complete 14-sublineage Jackknife Leave-One-Out (n = 13 per iteration)
    for s_out in df_m["sublineage"]:
        sub = df_m[df_m["sublineage"] != s_out]
        rho, p_sp, r, p_pe = correlations(sub, "n_countries", "retain_vs_2000")
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
            "P_Value_Primary": float(p_sp),
            "Statistic_Secondary": f"MWU_Global_vs_NonGlobal_P={p_mwu:.6f}",
            "P_Value_Secondary": float(p_mwu),
            "Significant_0.05": "YES" if (p_sp < 0.05 and p_mwu < 0.05) else "NO",
        })

    df_stats = pd.DataFrame(stat_rows)
    # BH-FDR separately for each correlation analysis set and method,
    # and for the nine main-cohort MWU contrasts. LOO is descriptive.
    df_stats["Q_Value_Primary_BH"] = np.nan
    df_stats["Q_Value_Secondary_BH"] = np.nan
    df_stats["Significant_Primary_FDR_0.05"] = "NA"
    for (_, kind), group in df_stats.groupby(["Analysis_Set", "Test_Type"]):
        if kind not in ("Correlation", "MannWhitneyU_TwoSided"):
            continue
        q = bh_adjust(group["P_Value_Primary"])
        df_stats.loc[group.index, "Q_Value_Primary_BH"] = q
        df_stats.loc[group.index, "Significant_Primary_FDR_0.05"] = np.where(q < 0.05, "YES", "NO")
        if kind == "Correlation":
            df_stats.loc[group.index, "Q_Value_Secondary_BH"] = bh_adjust(group["P_Value_Secondary"])
    df_m.to_csv(TABLES_DIR / "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv", sep="\t", index=False)
    df_stats.to_csv(TABLES_DIR / "Table_S2b_CrossFigure_Statistical_Tests.tsv", sep="\t", index=False)

    # Print key headline statistics for quick verification
    rho_main, p_main = stats.spearmanr(df_m["n_countries"], df_m["retain_vs_2000"])
    _, p_mwu_main = stats.mannwhitneyu(g_df["retain_vs_2000"], ng_df["retain_vs_2000"], alternative="two-sided")
    return df_m, df_stats


def main():
    run_cross_figure_statistics()


if __name__ == "__main__":
    main()
