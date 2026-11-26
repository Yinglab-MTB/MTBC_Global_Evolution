# -*- coding: utf-8 -*-
"""
Fig. 4 (Nature Editorial V2): Asymmetric antibiotic-era demographic resilience.
  a. Grand Normalized Overlay BSP Trajectory (1750-2025, Ne/Ne_peak) with shaded historical eras
  b. Compact 2x7 raw log-Ne sparkline grid with 95% HPD ribbons across all 14 sublineages
  c. Bubble scatter with 95% regression confidence band (n_countries vs 21st-century Ne retention)
  d. Editorial Raincloud Plot (half-violin + inner box + beeswarm): 21st-century Ne retention (P = 0.0020)
  e. Editorial Raincloud Plot contrast: overall MDR-TB prevalence (P = 0.364, n.s.; L2.1 spotlight)
"""

import os
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, CLEAN_DIR, PAL_GROUP,
    set_nature_v2_style, save_v2, draw_raincloud,
)


def load_bsp(sn):
    short = sn.replace("lineage", "L")
    if sn == "lineage4.1":
        p = os.path.join(CLEAN_DIR, "L4.1_clean499.BSP.txt")
    else:
        p = os.path.join(BASE_DIR, "Fig. 4", f"{short}.BSP.txt")
    lines = [l.strip() for l in open(p, encoding="utf-8", errors="ignore") if l.strip()]
    h_idx = next(i for i, l in enumerate(lines) if l.startswith("Time"))
    header = lines[h_idx].split("\t")
    rows = [[float(x) if x not in ("", "NaN", "nan") else np.nan for x in l.split("\t")] for l in lines[h_idx + 1:]]
    return pd.DataFrame(rows, columns=header).dropna(subset=["Median"]).sort_values("Time").reset_index(drop=True)


def main():
    print("[Fig. 4 V2] Generating Overlay BSP + Sparklines + CI Scatter + Raincloud Plots...")
    set_nature_v2_style()

    df_m = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv"), sep="\t")
    meta_idx = df_m.set_index("sublineage")

    fig = plt.figure(figsize=(16.5, 11.4))
    outer = gridspec.GridSpec(2, 1, height_ratios=[1.25, 1.0], hspace=0.30)

    # Top row: Left = Grand Normalized Overlay (a), Right = 2x7 Sparkline Grid (b)
    gs_top = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[1.15, 1.65], wspace=0.22)

    # --- Panel a: Grand Normalized Overlay BSP (1750-2025) ---
    ax_a = fig.add_subplot(gs_top[0])
    # Shaded historical eras
    ax_a.axvspan(1750, 1950, color="#FAF8F4", zorder=0)
    ax_a.axvspan(1950, 2000, color="#F3F6F8", zorder=0)
    ax_a.axvspan(2000, 2025, color="#FBECE8", alpha=0.45, zorder=0)
    ax_a.axvline(1950, color="#525B66", linestyle="--", linewidth=0.9, zorder=1)
    ax_a.axvline(2000, color="#C44536", linestyle=":", linewidth=0.9, zorder=1)

    ax_a.text(1850, 1.07, "Pre-antibiotic expansion\n(globalization & urbanization)", ha="center", va="top", fontsize=7.3, color="#6B7280", style="italic")
    ax_a.text(1975, 1.07, "Chemotherapy\nrollout", ha="center", va="top", fontsize=7.1, color="#525B66", style="italic")
    ax_a.text(2012.5, 1.07, "21st c.", ha="center", va="top", fontsize=7.1, color="#9B2C2C", fontweight="bold")

    # Plot Intermediate lightly first, then Rooted, then Global bold on top!
    order_draw = (
        [s for s in df_m["sublineage"] if meta_idx.loc[s, "ecological_group"] == "Intermediate"] +
        [s for s in df_m["sublineage"] if meta_idx.loc[s, "ecological_group"] == "Rooted"] +
        [s for s in df_m["sublineage"] if meta_idx.loc[s, "ecological_group"] == "Global"]
    )
    for sn in order_draw:
        d = load_bsp(sn)
        d = d[d["Time"] >= 1740].copy()
        grp = meta_idx.loc[sn, "ecological_group"]
        pk = meta_idx.loc[sn, "BSP_Ne_peak"]
        y_norm = d["Median"].values / pk
        col = PAL_GROUP[grp]
        lw = 2.1 if grp == "Global" else (2.0 if grp == "Rooted" else 1.0)
        alpha = 0.92 if grp in ("Global", "Rooted") else 0.35
        ls = "--" if grp == "Rooted" else "-"
        ax_a.plot(d["Time"], y_norm, color=col, linewidth=lw, alpha=alpha, linestyle=ls, zorder=4 if grp != "Intermediate" else 2)

    # Direct right-edge end labels for key Global & Rooted clades
    end_labels = [
        ("lineage4.8", 0.395, "L4.8"), ("lineage2.2", 0.355, "L2.2"), ("lineage4.1", 0.295, "L4.1"),
        ("lineage4.4", 0.265, "L4.4"), ("lineage4.3", 0.175, "L4.3"),
        ("lineage4.5", 0.095, "L4.5 (Rooted)"), ("lineage2.1", 0.042, "L2.1 (Rooted)"),
    ]
    for sn, y_e, lbl in end_labels:
        grp = meta_idx.loc[sn, "ecological_group"]
        ax_a.text(2027, y_e, lbl, va="center", fontsize=7.2, fontweight="bold", color=PAL_GROUP[grp])

    ax_a.set_xlim(1750, 2065)
    ax_a.set_ylim(0, 1.12)
    ax_a.set_xticks([1750, 1800, 1850, 1900, 1950, 2000, 2025])
    ax_a.set_xlabel("Calendar year (CE)")
    ax_a.set_ylabel("Normalized effective population size ($N_e / N_{e,\\mathrm{peak}}$)")
    ax_a.set_title("a   Normalized demographic trajectories across historical eras", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    handles_a = [
        plt.Line2D([0], [0], color=PAL_GROUP["Global"], lw=2.2, label="Global generalists (n = 5)"),
        plt.Line2D([0], [0], color=PAL_GROUP["Rooted"], lw=2.0, ls="--", label="Rooted specialists (n = 2)"),
        plt.Line2D([0], [0], color=PAL_GROUP["Intermediate"], lw=1.2, alpha=0.5, label="Intermediate (n = 7)"),
    ]
    ax_a.legend(handles=handles_a, loc="upper left", bbox_to_anchor=(0.01, 0.86), frameon=True,
                facecolor="white", edgecolor="#D6D3CD", fontsize=7.5)

    # --- Panel b: Compact 2x7 Individual Sparkline BSP Strip ---
    gs_sp = gridspec.GridSpecFromSubplotSpec(2, 7, subplot_spec=gs_top[1], hspace=0.45, wspace=0.22)
    row1 = ["lineage2.2", "lineage4.1", "lineage4.3", "lineage4.4", "lineage4.8", "lineage2.1", "lineage4.5"]
    row2 = ["lineage1.1", "lineage1.2", "lineage3.1", "lineage4.2", "lineage4.6", "lineage4.7", "lineage4.9"]
    for idx, sn in enumerate(row1 + row2):
        ax = fig.add_subplot(gs_sp[idx // 7, idx % 7])
        d = load_bsp(sn)
        meta = meta_idx.loc[sn]
        grp = meta["ecological_group"]
        col = PAL_GROUP[grp]
        ax.axvspan(1950, 2025, color="#FBECE8" if grp == "Global" else "#EEF2F5", alpha=0.6, zorder=0)
        ax.axvline(1950, color="#6B7280", ls="--", lw=0.6, zorder=1)
        ax.fill_between(d["Time"], np.clip(d["Lower"], 10, None), np.clip(d["Upper"], 10, None),
                        color=col, alpha=0.20, lw=0, zorder=2)
        ax.plot(d["Time"], np.clip(d["Median"], 10, None), color=col, lw=1.35, zorder=3)
        ax.set_yscale("log")
        ax.set_xlim(1450, 2028)
        ax.set_ylim(50, 1.8e6)
        ax.set_xticks([1500, 1950])
        ax.set_xticklabels(["1500", "1950"], fontsize=6.4)
        if idx % 7 == 0:
            ax.set_ylabel("log $N_e$", fontsize=7.5)
            ax.tick_params(axis="y", labelsize=6.3)
        else:
            ax.set_yticklabels([])
        short = sn.replace("lineage", "L")
        ax.set_title(f"{short}", loc="left", fontsize=8.0, fontweight="bold", color=col, pad=2)
        ax.text(0.06, 0.08, f"{meta['retain_vs_2000']*100:.0f}%", transform=ax.transAxes, fontsize=6.6,
                fontweight="bold", color=col, bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="#D6D3CD", alpha=0.9))

    fig.text(0.425, 0.905, "b   Individual sublineage Bayesian Skyline Plots (log $N_e$, 95% HPD; badge: 2000→2025 retention)",
             fontsize=10.3, fontweight="bold", color="#1E2229")

    # ==========================================================================
    # BOTTOM ROW: c (Bubble Scatter + 95% CI Band), d & e (Raincloud Plots)
    # ==========================================================================
    gs_bot = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=outer[1], width_ratios=[1.35, 1.0, 1.0], wspace=0.28)

    # --- Panel c: Bubble Scatter + 95% Regression CI Band ---
    ax_c = fig.add_subplot(gs_bot[0])
    x_v = df_m["n_countries"].values
    y_v = df_m["retain_vs_2000"].values * 100.0
    slope, inter, r_val, p_pe, std_err = stats.linregress(x_v, y_v)
    rho, p_sp = stats.spearmanr(x_v, y_v)
    xs = np.linspace(12, 92, 100)
    ys_hat = slope * xs + inter
    # Exact 95% t-confidence band around OLS fit
    n_pts = len(x_v)
    resid = y_v - (slope * x_v + inter)
    s_err = np.sqrt(np.sum(resid**2) / (n_pts - 2))
    t_crit = stats.t.ppf(0.975, n_pts - 2)
    ci_band = t_crit * s_err * np.sqrt(1.0 / n_pts + (xs - x_v.mean())**2 / np.sum((x_v - x_v.mean())**2))

    ax_c.fill_between(xs, ys_hat - ci_band, ys_hat + ci_band, color="#CBD5E0", alpha=0.38, zorder=1)
    ax_c.plot(xs, ys_hat, color="#2B2D42", lw=1.4, ls="--", zorder=2)

    for grp in ["Global", "Intermediate", "Rooted"]:
        sub = df_m[df_m["ecological_group"] == grp]
        sz = 42 + 120 * (np.log10(sub["BSP_Ne_peak"]) - 4.0)
        sz = np.clip(sz, 48, 210)
        ax_c.scatter(sub["n_countries"], sub["retain_vs_2000"] * 100, s=sz, color=PAL_GROUP[grp],
                     alpha=0.90, edgecolor="white", linewidth=0.9, label=grp, zorder=4)

    offsets = {
        "L4.1": (-7.2, 4.2), "L4.8": (1.5, 3.8), "L2.2": (-9.8, -6.2), "L4.3": (2.0, -5.2),
        "L4.4": (2.2, 1.5), "L4.2": (2.2, 1.0), "L3.1": (2.4, 0.5), "L1.2": (2.2, -5.0),
        "L4.7": (2.2, 1.0), "L1.1": (2.4, 1.0), "L4.6": (2.2, 0.5), "L4.5": (2.6, 1.8),
        "L2.1": (2.2, 1.5), "L4.9": (2.2, 1.0),
    }
    for _, r in df_m.iterrows():
        lbl = r["sublineage"].replace("lineage", "L")
        dx, dy = offsets.get(lbl, (2.0, 1.2))
        ax_c.text(r["n_countries"] + dx, r["retain_vs_2000"] * 100 + dy, lbl, fontsize=7.8, color="#1E2229")

    ax_c.set_xlim(9, 96)
    ax_c.set_ylim(-2, 128)
    ax_c.set_xlabel("Geographic breadth (number of countries)")
    ax_c.set_ylabel("21st-century $N_e$ retention ($N_{e,2025}/N_{e,2000}$, %)")
    ax_c.set_title("c   Geographic breadth predicts modern resilience", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_c.text(0.03, 0.96, f"Spearman $\\rho$ = {rho:+.3f}, P = {p_sp:.4f}\nPearson r = {r_val:+.3f}, P = {p_pe:.4f}\n(dot area scaled to $\\log_{{10}}$ peak $N_e$)",
              transform=ax_c.transAxes, va="top", fontsize=7.8,
              bbox=dict(boxstyle="round,pad=0.28", facecolor="#FAF8F4", edgecolor="#D6D3CD"))
    ax_c.legend(frameon=False, loc="lower right", fontsize=8.2)
    ax_c.grid(linestyle=":", alpha=0.3)

    # --- Panel d: Raincloud Plot (21st-century Ne retention) ---
    ax_d = fig.add_subplot(gs_bot[1])
    g_sub = df_m[df_m["ecological_group"] == "Global"]
    ng_sub = df_m[df_m["ecological_group"] != "Global"]
    g_ret = g_sub["retain_vs_2000"].values * 100
    ng_ret = ng_sub["retain_vs_2000"].values * 100
    _, p_mwu = stats.mannwhitneyu(g_ret, ng_ret, alternative="two-sided")

    draw_raincloud(
        ax_d,
        [g_ret, ng_ret],
        positions=[0.0, 1.0],
        colors=[PAL_GROUP["Global"], PAL_GROUP["Intermediate"]],
        dot_color_lists=[
            [PAL_GROUP["Global"]] * len(g_ret),
            [PAL_GROUP[g] for g in ng_sub["ecological_group"]],
        ],
    )
    # Significance bracket
    ax_d.plot([0, 0, 1, 1], [109, 113, 113, 109], color="#2B2D42", lw=0.9)
    ax_d.text(0.5, 114.5, f"Mann–Whitney P = {p_mwu:.4f}\nMedian: 98.0% vs 40.2%",
              ha="center", va="bottom", fontsize=7.9, fontweight="bold", color="#9B2C2C")
    ax_d.set_xticks([0, 1])
    ax_d.set_xticklabels(["Global\n(n = 5)", "Non-Global\n(n = 9)"], fontsize=8.5, fontweight="bold")
    ax_d.set_xlim(-0.48, 1.52)
    ax_d.set_ylim(-4, 130)
    ax_d.set_ylabel("21st-century $N_e$ retention (%)")
    ax_d.set_title("d   Modern demographic resilience", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_d.grid(axis="y", linestyle=":", alpha=0.3)

    # --- Panel e: Raincloud Plot Contrast (Overall MDR-TB Rate) ---
    ax_e = fig.add_subplot(gs_bot[2])
    g_mdr = g_sub["MDR_rate"].values * 100
    ng_mdr = ng_sub["MDR_rate"].values * 100
    _, p_mdr = stats.mannwhitneyu(g_mdr, ng_mdr, alternative="two-sided")

    draw_raincloud(
        ax_e,
        [g_mdr, ng_mdr],
        positions=[0.0, 1.0],
        colors=[PAL_GROUP["Global"], PAL_GROUP["Intermediate"]],
        dot_color_lists=[
            [PAL_GROUP["Global"]] * len(g_mdr),
            [PAL_GROUP[g] for g in ng_sub["ecological_group"]],
        ],
    )
    ax_e.plot([0, 0, 1, 1], [47.5, 49.5, 49.5, 47.5], color="#525B66", lw=0.9)
    ax_e.text(0.5, 50.5, f"P = {p_mdr:.3f} (n.s.)\nMedian: 22.2% vs 18.3%",
              ha="center", va="bottom", fontsize=7.9, color="#525B66")
    ax_e.annotate("L2.1: 44.3% MDR,\nyet rooted & collapsed", xy=(0.82, 44.3), xytext=(0.10, 39.2),
                  fontsize=7.3, fontweight="bold", color="#525B66",
                  bbox=dict(boxstyle="round,pad=0.2", facecolor="#FAF8F4", edgecolor="#CBD5E0"),
                  arrowprops=dict(arrowstyle="-|>", color="#525B66", lw=0.9))
    ax_e.set_xticks([0, 1])
    ax_e.set_xticklabels(["Global\n(n = 5)", "Non-Global\n(n = 9)"], fontsize=8.5, fontweight="bold")
    ax_e.set_xlim(-0.48, 1.52)
    ax_e.set_ylim(6, 58)
    ax_e.set_ylabel("Overall MDR-TB prevalence (%)")
    ax_e.set_title("e   Total MDR burden alone ≠ global success", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_e.grid(axis="y", linestyle=":", alpha=0.3)

    save_v2(fig, "Fig4_Demographic_Resilience_BSP_V2")


if __name__ == "__main__":
    main()
