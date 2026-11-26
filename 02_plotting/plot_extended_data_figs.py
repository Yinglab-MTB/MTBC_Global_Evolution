# -*- coding: utf-8 -*-
"""
Extended Data Figs. 1-4 (Nature Editorial V2):
  - Ext. Data Fig. 1 V2: Sankey-style curved-ribbon Cohort Funnel + editorial DR & Lineage strips
  - Ext. Data Fig. 2 V2: Shaded Ecological Niche Quadrant + Grouped Lollipop of 30 L3 sublineages
  - Ext. Data Fig. 3 V2: Lollipop BEAST QC Audit + Leave-L4.1-out CI Sensitivity + Clade Artifact Strip
  - Ext. Data Fig. 4 V2: Bubble-Heatmap (14x10 DR) + Paired Slopegraph + Diverging L2.2/L2.1 Butterfly + L4 Lollipop
"""

import os
import textwrap
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.path import Path
from matplotlib.patches import PathPatch, FancyBboxPatch, Patch

import sys
from pathlib import Path as _FilePath
sys.path.insert(0, str(_FilePath(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, PAL_GROUP, PAL_LINEAGE, CMAP_BORDEAUX,
    set_nature_v2_style, save_v2,
)


def draw_sankey_band(ax, y_top, w_top, y_bot, w_bot, color, alpha=0.22):
    x0_t, x1_t = 5.0 - w_top / 2.0, 5.0 + w_top / 2.0
    x0_b, x1_b = 5.0 - w_bot / 2.0, 5.0 + w_bot / 2.0
    ym = 0.5 * (y_top + y_bot)
    verts = [
        (x0_t, y_top), (x1_t, y_top),
        (x1_t, ym), (x1_b, ym), (x1_b, y_bot),
        (x0_b, y_bot),
        (x0_b, ym), (x0_t, ym), (x0_t, y_top),
    ]
    codes = [
        Path.MOVETO, Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
        Path.LINETO,
        Path.CURVE4, Path.CURVE4, Path.CURVE4,
    ]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=1))


def make_ext_fig1_v2():
    print("[Ext. Data Fig. 1 V2] Generating Sankey Cohort Funnel & Editorial Taxonomic/DR Strips...")
    set_nature_style_fn = set_nature_v2_style
    set_nature_style_fn()

    df_f = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S1_Cohort_Filtering_Funnel.tsv"), sep="\t")
    df_b = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S1b_MainLineage_and_DR_Breakdown.tsv"), sep="\t")
    df_dr = df_b[df_b["Category"] == "Drug_Resistance_Type"].sort_values("Count", ascending=True)
    df_lin = df_b[df_b["Category"] == "Main_Lineage"].sort_values("Count", ascending=True)

    fig = plt.figure(figsize=(16.0, 8.6))
    gs = gridspec.GridSpec(2, 2, width_ratios=[1.25, 1.0], hspace=0.34, wspace=0.28)

    # Panel a: Sankey-Style Curved Ribbon Cohort Funnel
    ax_a = fig.add_subplot(gs[:, 0])
    ax_a.set_xlim(0, 10)
    ax_a.set_ylim(0, 11.2)
    ax_a.axis("off")
    ax_a.set_title("a   Four-stage hierarchical cohort filtering funnel", loc="left", fontweight="bold", fontsize=10.8, pad=6)

    pal_f = ["#1B365D", "#2A6F82", "#3E8E7E", "#D99B26", "#C44536"]
    bg_f  = ["#EEF3F8", "#EDF6F8", "#EFF7F4", "#FDF8EC", "#FDF2EF"]
    widths = [9.2, 8.7, 8.2, 7.8, 7.4]
    ys = [9.2, 7.0, 4.8, 2.6, 0.4]
    h_box = 1.48

    for i in range(len(df_f) - 1):
        draw_sankey_band(ax_a, ys[i], widths[i] * 0.72, ys[i + 1] + h_box, widths[i + 1] * 0.72, pal_f[i + 1], alpha=0.22)

    for i, (_, r) in enumerate(df_f.iterrows()):
        w = widths[i]
        x = (10.0 - w) / 2.0
        y = ys[i]
        ax_a.add_patch(FancyBboxPatch((x, y), w, h_box, boxstyle="round,pad=0.12,rounding_size=0.22",
                                      facecolor=bg_f[i], edgecolor=pal_f[i], lw=1.4, zorder=2))
        ax_a.add_patch(FancyBboxPatch((x - 0.05, y + 0.06), 0.16, h_box - 0.12, boxstyle="round,pad=0.02,rounding_size=0.08",
                                      facecolor=pal_f[i], edgecolor="none", zorder=3))
        ax_a.text(x + 0.30, y + h_box - 0.28, f"{r['Stage']}   (n = {int(r['N_Genomes']):,})",
                  ha="left", va="top", fontsize=8.3, fontweight="bold", color=pal_f[i], zorder=4)
        scope = "\n".join(textwrap.wrap(str(r["Scope_Details"]), width=68))
        ax_a.text(x + 0.30, y + 0.44, scope, ha="left", va="center", fontsize=6.9, color="#2B2D42", linespacing=1.22, zorder=4)

    # Panel b: Editorial Lollipop WHO Drug-Resistance Breakdown
    ax_b = fig.add_subplot(gs[0, 1])
    dr_pal = {
        "Sensitive": "#2A6F82", "Other": "#A89F91", "HR-TB": "#D99B26",
        "RR-TB": "#E59866", "MDR-TB": "#C44536", "Pre-XDR-TB": "#7A1C35", "XDR-TB": "#380920",
    }
    y_d = np.arange(len(df_dr))
    cols_d = [dr_pal.get(k, "#6B7280") for k in df_dr["Label"]]
    ax_b.hlines(y_d, 0, df_dr["Count"], color=cols_d, lw=4.5, alpha=0.82, zorder=2)
    ax_b.scatter(df_dr["Count"], y_d, s=75, color=cols_d, edgecolor="white", lw=0.8, zorder=3)
    ax_b.set_yticks(y_d)
    ax_b.set_yticklabels(df_dr["Label"], fontsize=8.4, fontweight="bold")
    ax_b.set_xlim(0, 105000)
    for yi, (_, r) in zip(y_d, df_dr.iterrows()):
        ax_b.text(r["Count"] + 1800, yi, f"{int(r['Count']):,} ({r['Percentage']:.1f}%)", va="center", fontsize=7.6, color="#2B2D42")
    ax_b.set_xlabel("Number of genomes (total n = 162,869)")
    ax_b.set_title("b   Genotypic drug-resistance categories", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_b.grid(axis="x", linestyle=":", alpha=0.35)

    # Panel c: Editorial Lollipop Complete Lineage Counts (log scale)
    ax_c = fig.add_subplot(gs[1, 1])
    y_l = np.arange(len(df_lin))
    cols_l = [PAL_LINEAGE.get(k, "#6B7280") for k in df_lin["Label"]]
    ax_c.hlines(y_l, 1, df_lin["Count"], color=cols_l, lw=3.2, alpha=0.82, zorder=2)
    ax_c.scatter(df_lin["Count"], y_l, s=56, color=cols_l, edgecolor="white", lw=0.7, zorder=3)
    ax_c.set_xscale("log")
    ax_c.set_xlim(1, 7.5e5)
    ax_c.set_yticks(y_l)
    ax_c.set_yticklabels([k.replace("lineage", "L") for k in df_lin["Label"]], fontsize=8.0, fontweight="bold")
    for yi, (_, r) in zip(y_l, df_lin.iterrows()):
        ax_c.text(r["Count"] * 1.25, yi, f"{int(r['Count']):,} ({r['Percentage']:.2f}%)", va="center", fontsize=7.2, color="#2B2D42")
    ax_c.set_xlabel("Number of genomes (log scale)")
    ax_c.set_title("c   Complete MTBC lineage composition", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_c.grid(axis="x", linestyle=":", alpha=0.35)

    save_v2(fig, "Extended_Data_Fig1_Cohort_Funnel_V2")


def make_ext_fig2_v2():
    print("[Ext. Data Fig. 2 V2] Generating Ecological Niche Quadrant + L3 Lollipop Ranking...")
    set_nature_v2_style()

    df_s2 = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv"), sep="\t")
    df_s3 = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2c_ThirdLevel_Sublineage_Geography.tsv"), sep="\t")

    fig = plt.figure(figsize=(16.4, 8.2))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 1.15], wspace=0.52)

    # Panel a: Shaded Ecological Niche Quadrant
    ax_a = fig.add_subplot(gs[0])
    ax_a.axvspan(62, 99, ymin=0.635, ymax=1.0, color="#FBECE8", alpha=0.55, zorder=0)
    ax_a.axvspan(8, 40, ymin=0.0, ymax=0.23, color="#EEF0F2", alpha=0.75, zorder=0)
    ax_a.axvline(62, color="#C44536", ls="--", lw=0.85, alpha=0.65)
    ax_a.axhline(2.75, color="#C44536", ls="--", lw=0.85, alpha=0.65)
    ax_a.axhline(2.00, color="#6B7280", ls="--", lw=0.85, alpha=0.65)

    for grp in ["Global", "Intermediate", "Rooted"]:
        sub = df_s2[df_s2["ecological_group"] == grp]
        sz = 50 + 220 * np.sqrt(sub["N_geo"] / df_s2["N_geo"].max())
        ax_a.scatter(sub["n_countries"], sub["shannon_H"], s=sz, color=PAL_GROUP[grp],
                     alpha=0.90, edgecolor="white", lw=1.0, label=f"{grp} (n={len(sub)})", zorder=4)

    offsets = {
        "L2.2": (-7.8, -0.09), "L4.3": (-9.5, -0.075), "L4.1": (1.8, 0.03), "L4.8": (-15.0, 0.03),
        "L4.4": (2.2, -0.015), "L4.9": (1.8, 0.03), "L4.2": (1.8, -0.03), "L4.7": (1.8, 0.03),
        "L4.6": (-8.0, 0.045), "L1.2": (1.8, 0.02), "L1.1": (1.8, -0.04), "L3.1": (1.8, 0.02),
        "L4.5": (1.8, 0.03), "L2.1": (1.8, 0.03),
    }
    for _, r in df_s2.iterrows():
        lbl = r["sublineage"].replace("lineage", "L")
        dx, dy = offsets.get(lbl, (1.8, 0.02))
        ax_a.text(r["n_countries"] + dx, r["shannon_H"] + dy, f"{lbl} ({r['top3_share']*100:.0f}%)",
                  fontsize=7.6, color="#1E2229", fontweight="bold" if r["ecological_group"] == "Global" else "normal")

    ax_a.text(64, 3.29, "Cosmopolitan Generalist Niche\n(≥65 countries, H > 2.75)", fontsize=7.6, color="#8C261D", fontweight="bold")
    ax_a.text(10, 1.65, "Rooted Specialist Niche\n(H < 2.0, Top-3 ≥ 75%)", fontsize=7.6, color="#4A5568", fontweight="bold")
    ax_a.set_xlim(8, 99)
    ax_a.set_ylim(1.60, 3.38)
    ax_a.set_xlabel("Geographic breadth (number of countries)")
    ax_a.set_ylabel("Shannon geographic entropy (H)")
    ax_a.set_title("a   Quantitative ecological niche classification (14 L2 sublineages)", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_a.legend(loc="lower right", frameon=True, facecolor="#FAF8F4", edgecolor="#D6D3CD", fontsize=8.0)
    ax_a.grid(linestyle=":", alpha=0.3)

    # Panel b: Editorial Lollipop Ranking of all 30 Third-Level Sublineages
    ax_b = fig.add_subplot(gs[1])
    df_s3_ord = df_s3.sort_values("n_countries", ascending=True).reset_index(drop=True)
    y = np.arange(len(df_s3_ord))
    cols = [PAL_GROUP[g] for g in df_s3_ord["parent_group"]]
    ax_b.hlines(y, 0, df_s3_ord["n_countries"], color=cols, lw=2.4, alpha=0.75, zorder=2)
    ax_b.scatter(df_s3_ord["n_countries"], y, s=38, color=cols, edgecolor="white", lw=0.6, zorder=3)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([f"{s.replace('lineage','L')}  [{t}]" for s, t in zip(df_s3_ord["sublineage_L3"], df_s3_ord["top1_country"])], fontsize=7.1)
    for yi, (_, r) in zip(y, df_s3_ord.iterrows()):
        h_val = abs(float(r["shannon_H"]))
        ax_b.text(r["n_countries"] + 1.8, yi, f"{int(r['n_countries'])} (H={h_val:.2f})", va="center", fontsize=6.7, color="#2B2D42")
    ax_b.set_xlim(0, 110)
    ax_b.set_xlabel("Countries colonized (label: dominant country & share)")
    ax_b.set_title("b   Ranked breadth across all 30 third-level sublineages", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_b.grid(axis="x", linestyle=":", alpha=0.3)

    save_v2(fig, "Extended_Data_Fig2_Spatial_Metrics_L2_L3_V2")


def make_ext_fig3_v2():
    print("[Ext. Data Fig. 3 V2] Generating Editorial QC, Leave-L4.1-out Sensitivity & Clade Artifact Strip...")
    set_nature_v2_style()

    df_qc = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S4_BEAST_Tree_QC_and_Outlier_Audit.tsv"), sep="\t")
    df_s2 = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv"), sep="\t")
    df_art = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S3b_Flagged_Clade_Hitchhiking_Candidates.tsv"), sep="\t")

    fig = plt.figure(figsize=(16.2, 5.4))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1.12, 1.0, 1.15], wspace=0.30)

    # Panel a: Lollipop BEAST Tree QC Ratio (H1/H2) — Post-QC after removing L4.1 outlier 2405024009
    ax_a = fig.add_subplot(gs[0])
    df_qcs = df_qc.sort_values("H1_H2_ratio", ascending=True).reset_index(drop=True)
    y = np.arange(len(df_qcs))
    cols = ["#C44536" if r >= 1.8 else ("#D99B26" if r >= 1.35 else "#2A6F82") for r in df_qcs["H1_H2_ratio"]]
    ax_a.axvspan(1.80, 2.25, color="#FBECE8", alpha=0.5, zorder=0)
    ax_a.axvline(1.35, color="#D99B26", ls="--", lw=0.85, label="Basal check (1.35)")
    ax_a.axvline(1.80, color="#C44536", ls="--", lw=0.95, label="Outlier threshold (1.80)")
    ax_a.hlines(y, 1.0, df_qcs["H1_H2_ratio"], color=cols, lw=2.8, alpha=0.82, zorder=2)
    ax_a.scatter(df_qcs["H1_H2_ratio"], y, s=54, color=cols, edgecolor="white", lw=0.7, zorder=3)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels([s.replace("lineage", "L") for s in df_qcs["sublineage"]], fontsize=8.0, fontweight="bold")
    for yi, (_, r) in zip(y, df_qcs.iterrows()):
        is_l41 = (r["sublineage"] == "lineage4.1")
        extra = "  (post-QC; 2405024009 removed)" if is_l41 else ""
        bbox_kw = dict(facecolor="white", edgecolor="none", alpha=0.92, pad=1.0) if is_l41 else None
        ax_a.text(r["H1_H2_ratio"] + 0.03, yi, f"{r['H1_H2_ratio']:.2f}{extra}", va="center",
                  fontsize=7.0, fontweight="bold" if is_l41 else "normal",
                  color="#2A6F82" if is_l41 else "#2B2D42", bbox=bbox_kw, zorder=5)
    ax_a.set_xlim(0.95, 2.22)
    ax_a.set_xlabel("Root / second-deepest node height ratio ($H_1/H_2$)")
    ax_a.set_title("a   Post-QC BEAST tree balance (all $H_1/H_2$ < 1.80)", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_a.legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#D6D3CD", fontsize=7.3)

    # Panel b: Jackknife Leave-One-Out Sensitivity (Original Scatter + 95% CI + 14 LOO fits)
    ax_b = fig.add_subplot(gs[1])
    xv_all = df_s2["n_countries"].values
    yv_all = df_s2["retain_vs_2000"].values * 100
    slope_all, inter_all, _, _, _ = stats.linregress(xv_all, yv_all)
    xs = np.linspace(12, 92, 100)
    ys_hat = slope_all * xs + inter_all
    s_err = np.sqrt(np.sum((yv_all - (slope_all * xv_all + inter_all)) ** 2) / (len(xv_all) - 2))
    ci = stats.t.ppf(0.975, len(xv_all) - 2) * s_err * np.sqrt(
        1 / len(xv_all) + (xs - xv_all.mean()) ** 2 / np.sum((xv_all - xv_all.mean()) ** 2)
    )

    ax_b.fill_between(xs, ys_hat - ci, ys_hat + ci, color="#CBD5E0", alpha=0.38, zorder=1)
    # Overlay all 14 leave-one-out regression lines tightly inside the 95% CI band
    rhos, ps_sp, rs_pe, ps_pe, ps_mwu = [], [], [], [], []
    for idx_out, s_out in enumerate(df_s2["sublineage"]):
        sub = df_s2[df_s2["sublineage"] != s_out]
        xi, yi = sub["n_countries"].values, sub["retain_vs_2000"].values * 100
        sl_i, in_i, r_i, p_pe_i, _ = stats.linregress(xi, yi)
        rho_i, p_sp_i = stats.spearmanr(xi, yi)
        g_i = sub[sub["ecological_group"] == "Global"]["retain_vs_2000"]
        ng_i = sub[sub["ecological_group"] != "Global"]["retain_vs_2000"]
        _, p_mwu_i = stats.mannwhitneyu(g_i, ng_i, alternative="two-sided")
        rhos.append(rho_i); ps_sp.append(p_sp_i); rs_pe.append(r_i); ps_pe.append(p_pe_i); ps_mwu.append(p_mwu_i)
        ax_b.plot(xs, sl_i * xs + in_i, color="#64748B", lw=0.75, alpha=0.36, zorder=2,
                  label="14 LOO fits" if idx_out == 0 else None)

    ax_b.plot(xs, ys_hat, color="#2B2D42", ls="--", lw=1.3, zorder=3)
    for grp in ["Global", "Intermediate", "Rooted"]:
        s = df_s2[df_s2["ecological_group"] == grp]
        ax_b.scatter(s["n_countries"], s["retain_vs_2000"] * 100, s=66, color=PAL_GROUP[grp],
                     edgecolor="white", lw=0.8, label=grp, zorder=4)

    ax_b.set_xlim(10, 95)
    ax_b.set_ylim(-2, 132)
    ax_b.set_xlabel("Geographic breadth (number of countries)")
    ax_b.set_ylabel("21st-century $N_e$ retention (%)")
    ax_b.set_title("b   Jackknife leave-one-out sensitivity (14 sets)", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_b.text(0.03, 0.96,
              f"All 14 LOO sets (n = 13 each):\n"
              f"Spearman $\\rho$ = {min(rhos):+.3f} to {max(rhos):+.3f} (P ≤ {max(ps_sp):.3f})\n"
              f"Pearson r = {min(rs_pe):+.3f} to {max(rs_pe):+.3f} (P ≤ {max(ps_pe):.3f})\n"
              f"Global vs Non-Global MWU: P ≤ {max(ps_mwu):.4f}",
              transform=ax_b.transAxes, va="top", fontsize=7.2,
              bbox=dict(boxstyle="round,pad=0.26", facecolor="#FAF8F4", edgecolor="#D6D3CD"))
    ax_b.legend(loc="lower right", frameon=False, fontsize=7.4)
    ax_b.grid(linestyle=":", alpha=0.3)

    # Panel c: Editorial Cluster Strip of 18 Filtered Clade-Hitchhiking Artifacts
    ax_c = fig.add_subplot(gs[2])
    df_arts = df_art.sort_values(["max_delta_freq", "gene_name"]).reset_index(drop=True)
    y_c = np.arange(len(df_arts))
    c_map = ["#6D597A" if "ClusterA" in d else ("#D99B26" if "ClusterB" in d else "#2A6F82") for d in df_arts["Audit_Diagnosis"]]
    ax_c.hlines(y_c, 0, df_arts["max_delta_freq"], color=c_map, lw=1.8, alpha=0.82, zorder=2)
    ax_c.scatter(df_arts["max_delta_freq"], y_c, color=c_map, s=44, edgecolor="white", lw=0.6, zorder=3)
    for xv_c, cc in [(0.2469, "#6D597A"), (0.3227, "#D99B26"), (0.6141, "#2A6F82")]:
        ax_c.axvline(xv_c, color=cc, ls=":", lw=0.95)
    ax_c.set_yticks(y_c)
    ax_c.set_yticklabels([f"{g}{'*' if 'VERIFY' in d else ''}" for g, d in zip(df_arts["gene_name"], df_arts["Audit_Diagnosis"])],
                         fontsize=7.4, fontstyle="italic")
    ax_c.set_xlim(0, 0.72)
    ax_c.set_xlabel("Identical single-SNP $\\Delta f$ across nested L4 contrasts")
    ax_c.set_title("c   Filtered clonal clade-hitchhiking artifacts (n = 18)", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    handles = [
        Patch(facecolor="#6D597A", label="Cluster A: L4.7–4.9/H37Rv (Δf=0.247, n=13)"),
        Patch(facecolor="#D99B26", label="Cluster B: L4.6–4.9 clade (Δf=0.323, n=3)"),
        Patch(facecolor="#2A6F82", label="Cluster C: L4-dominant (Δf=0.614, n=2)"),
    ]
    ax_c.legend(handles=handles, loc="lower right", frameon=True, facecolor="#FAF8F4", edgecolor="#D6D3CD", fontsize=6.9)

    save_v2(fig, "Extended_Data_Fig3_QC_Sensitivity_Hitchhiking_V2")


def make_ext_fig4_v2():
    print("[Ext. Data Fig. 4 V2] Generating Bubble-Heatmap + Slopegraph + Diverging Butterfly + L4 Lollipop...")
    set_nature_v2_style()

    dr_dir = os.path.join(BASE_DIR, "Extended Data Fig. 4")
    df_bur = pd.read_csv(os.path.join(dr_dir, "sublineage_resistance_burden.tsv"), sep="\t")
    df_drg = pd.read_csv(os.path.join(dr_dir, "DR_mutation_global_enriched.tsv"), sep="\t")
    df_wl = pd.read_csv(os.path.join(dr_dir, "DR_mutation_within_lineage_target_enriched.tsv"), sep="\t")
    df_s2 = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv"), sep="\t").set_index("sublineage")

    fig = plt.figure(figsize=(16.4, 11.2))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.12, 1.0], width_ratios=[1.25, 1.0], hspace=0.34, wspace=0.28)

    # --- Panel a: Bubble-Matrix of 14 Sublineages x 10 DR Phenotypes ---
    ax_a = fig.add_subplot(gs[0, 0])
    df_bur = df_bur.sort_values("MDR_rate", ascending=False).reset_index(drop=True)
    rate_cols = ["any_DR_rate", "INH_R_rate", "RIF_R_rate", "MDR_rate", "pre_XDR_rate",
                 "FQ_R_rate", "SLID_R_rate", "EMB_R_rate", "PZA_R_rate", "STR_R_rate"]
    short_cols = ["Any DR", "INH", "RIF", "MDR", "pre-XDR", "FQ", "SLID", "EMB", "PZA", "STR"]

    for i in range(len(df_bur)):
        if i % 2 == 1:
            ax_a.axhspan(i - 0.5, i + 0.5, color="#FAF8F4", zorder=0)
        grp = df_s2.loc[df_bur.loc[i, "sublineage"], "ecological_group"]
        ax_a.scatter([-0.75], [i], s=75, marker="s", color=PAL_GROUP[grp], edgecolor="none", zorder=3)
        for j, rc in enumerate(rate_cols):
            val = float(df_bur.loc[i, rc]) * 100.0
            norm_v = np.clip(val / 62.0, 0.0, 1.0)
            sz = 18 + 230 * (norm_v ** 1.15)
            ax_a.scatter([j], [i], s=sz, color=CMAP_BORDEAUX(norm_v),
                         edgecolor="#333840" if val >= 25 else "#CBD5E0", lw=0.55 if val >= 25 else 0.3, zorder=4)
            if val >= 18.0:
                ax_a.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=6.2,
                          fontweight="bold", color="white" if val >= 32 else "#1E2229", zorder=5)

    ax_a.set_xlim(-1.1, len(short_cols) - 0.35)
    ax_a.set_ylim(len(df_bur) - 0.5, -0.5)
    ax_a.set_xticks(range(len(short_cols)))
    ax_a.set_xticklabels(short_cols, fontsize=8.0, rotation=30, ha="right")
    ax_a.set_yticks(range(len(df_bur)))
    ax_a.set_yticklabels([f"{r['sublineage'].replace('lineage','L')} (n={int(r['n_samples']):,})" for _, r in df_bur.iterrows()], fontsize=7.9)
    ax_a.spines["left"].set_visible(False)
    ax_a.tick_params(axis="y", length=0)
    sm = plt.cm.ScalarMappable(cmap=CMAP_BORDEAUX, norm=plt.Normalize(0, 62))
    cb = fig.colorbar(sm, ax=ax_a, shrink=0.75, pad=0.02, aspect=24)
    cb.set_label("Resistance prevalence (%)", fontsize=8.2)
    ax_a.set_title("a   Sublineage drug-resistance matrix across 10 phenotypes (n = 23,533)", loc="left", fontweight="bold", fontsize=10.3, pad=6)

    # --- Panel b: Editorial Paired Slopegraph (Rooted -> Global) for 6 Canonical Low-Cost Mutations ---
    ax_b = fig.add_subplot(gs[0, 1])
    df_drg = df_drg.sort_values("freq_global", ascending=False).reset_index(drop=True)
    mut_cols = ["#C44536", "#7A1C35", "#2A6F82", "#D99B26", "#6D597A", "#3E8E7E"]

    ax_b.axvline(0, color="#CBD5E0", lw=1.2, zorder=1)
    ax_b.axvline(1, color="#CBD5E0", lw=1.2, zorder=1)

    for idx, (_, r) in enumerate(df_drg.iterrows()):
        fr, fg = r["freq_rooted"] * 100.0, r["freq_global"] * 100.0
        c = mut_cols[idx % len(mut_cols)]
        ax_b.plot([0, 1], [fr, fg], color=c, lw=2.2, alpha=0.85, zorder=3)
        ax_b.scatter([0, 1], [fr, fg], s=64, color=c, edgecolor="white", lw=0.8, zorder=4)
        # Left label
        y_l = fr
        if r["gene_name"] == "embB":
            y_l = fr - 0.5
        elif r["gene_name"] == "rpsL":
            y_l = fr + 0.3
        ax_b.text(-0.04, y_l, f"{fr:.1f}%", ha="right", va="center", fontsize=7.4, color="#525B66")
        # Right label
        y_r = fg
        if r["gene_name"] == "inhA":
            y_r = fg + 0.7
        elif r["gene_name"] == "embB":
            y_r = fg - 0.7
        fdr_s = f"{r['FDR']:.1e}" if r["FDR"] < 1e-3 else f"{r['FDR']:.3f}"
        ax_b.text(1.04, y_r, f"{fg:.1f}%  {r['gene_name']} {r['change']} (OR {r['odds_ratio']:.2f})",
                  ha="left", va="center", fontsize=7.5, fontweight="bold", color=c)

    ax_b.set_xlim(-0.28, 2.08)
    ax_b.set_ylim(-1.0, 31.0)
    ax_b.set_xticks([0, 1])
    ax_b.set_xticklabels(["Rooted\n(n = 941)", "Global\n(n = 15,257)"], fontsize=8.8, fontweight="bold")
    ax_b.set_ylabel("Mutation prevalence (%)")
    ax_b.spines["bottom"].set_visible(False)
    ax_b.set_title("b   Slopegraph of 6 low-fitness-cost mutations (Rooted → Global)", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_b.grid(axis="y", linestyle=":", alpha=0.3)

    # --- Panel c: Diverging Butterfly Chart (L2.1 left vs L2.2 right) ---
    ax_c = fig.add_subplot(gs[1, 0])
    r22 = df_bur[df_bur["sublineage"] == "lineage2.2"].iloc[0]
    r21 = df_bur[df_bur["sublineage"] == "lineage2.1"].iloc[0]
    l2_wl = df_wl[df_wl["comparison"] == "L22_vs_L21"].set_index("mutation_id")

    rows_bf = [
        ("Any DR", r21["any_DR_rate"]*100, r22["any_DR_rate"]*100, ""),
        ("MDR-TB", r21["MDR_rate"]*100, r22["MDR_rate"]*100, "Nearly equal MDR"),
        ("pre-XDR", r21["pre_XDR_rate"]*100, r22["pre_XDR_rate"]*100, ""),
        ("Fluoroquinolone-R", r21["FQ_R_rate"]*100, r22["FQ_R_rate"]*100, ""),
        ("embB p.Met306Val", l2_wl.loc["embB|p.Met306Val|ethambutol", "freq_control"]*100,
         l2_wl.loc["embB|p.Met306Val|ethambutol", "freq_target"]*100, "OR = 3.10 (FDR 1.1e-3)"),
        ("rpsL p.Lys43Arg", l2_wl.loc["rpsL|p.Lys43Arg|streptomycin", "freq_control"]*100,
         l2_wl.loc["rpsL|p.Lys43Arg|streptomycin", "freq_target"]*100, "OR = 10.95 (FDR 9.3e-24)"),
    ][::-1]
    y_bf = np.arange(len(rows_bf))
    ax_c.axvline(0, color="#2B2D42", lw=1.0, zorder=3)
    for yi, (nm, v21, v22, note) in enumerate(rows_bf):
        ax_c.barh(yi, -v21, height=0.58, color=PAL_GROUP["Rooted"], edgecolor="white", lw=0.6, alpha=0.88, zorder=2)
        ax_c.barh(yi, v22, height=0.58, color=PAL_GROUP["Global"], edgecolor="white", lw=0.6, alpha=0.90, zorder=2)
        ax_c.text(-v21 - 1.5, yi, f"{v21:.1f}%", ha="right", va="center", fontsize=7.6, color="#4A5568", fontweight="bold")
        txt_r = f"{v22:.1f}%  {note}" if note else f"{v22:.1f}%"
        ax_c.text(v22 + 1.5, yi, txt_r, ha="left", va="center", fontsize=7.5,
                  color="#8C261D" if "OR" in note else "#2B2D42", fontweight="bold")

    ax_c.set_yticks(y_bf)
    ax_c.set_yticklabels([r[0] for r in rows_bf], fontsize=8.2, fontweight="bold")
    ax_c.set_xlim(-78, 92)
    ax_c.set_xticks([-60, -40, -20, 0, 20, 40, 60])
    ax_c.set_xticklabels(["60%", "40%", "20%", "0", "20%", "40%", "60%"])
    ax_c.set_xlabel("← L2.1 Rooted (n = 262)   |   L2.2 Global (n = 4,288) →")
    ax_c.set_title("c   Diverging sister-clade contrast within Lineage 2", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_c.grid(axis="x", linestyle=":", alpha=0.35)

    # --- Panel d: Editorial Lollipop of Top Within-L4 Epistatic DR Mutations ---
    ax_d = fig.add_subplot(gs[1, 1])
    l4_comps = ["L41_vs_otherL4", "L43_vs_otherL4", "L44_vs_otherL4", "L48_vs_otherL4"]
    df_top4 = pd.concat([df_wl[df_wl["comparison"] == c].sort_values("FDR").head(3) for c in l4_comps], ignore_index=True)
    df_top4["neg_log10_FDR"] = -np.log10(np.clip(df_top4["FDR"].values, 1e-65, 1.0))
    df_top4["short_comp"] = df_top4["comparison"].str.replace("_vs_otherL4", "").str.replace("L4", "L4.")
    df_top4["lbl"] = df_top4["short_comp"] + " · " + df_top4["gene_name"] + " " + df_top4["change"]

    comp_pal = {"L4.1": "#C44536", "L4.3": "#2A6F82", "L4.4": "#3E8E7E", "L4.8": "#6D597A"}
    y_d = np.arange(len(df_top4))[::-1]
    c_d = [comp_pal[c] for c in df_top4["short_comp"]]
    ax_d.hlines(y_d, 0, df_top4["neg_log10_FDR"], color=c_d, lw=2.2, alpha=0.82, zorder=2)
    ax_d.scatter(df_top4["neg_log10_FDR"], y_d, s=35 + df_top4["freq_target"] * 1250,
                 color=c_d, edgecolor="white", lw=0.7, zorder=3)
    for yi, (_, r) in zip(y_d, df_top4.iterrows()):
        or_txt = "OR=∞" if np.isinf(r["odds_ratio"]) else f"OR={r['odds_ratio']:.1f}"
        ax_d.text(r["neg_log10_FDR"] + 1.4, yi, f"{r['freq_target']*100:.1f}% ({or_txt})", va="center", fontsize=7.1, color="#2B2D42")
    ax_d.set_yticks(y_d)
    ax_d.set_yticklabels(df_top4["lbl"], fontsize=7.6)
    ax_d.set_xlim(0, 76)
    ax_d.set_xlabel("$-\\log_{10}(\\mathrm{FDR})$ (dot size scaled to target prevalence)")
    ax_d.set_title("d   Within-L4 sublineage-enriched DR spectra (170 total)", loc="left", fontweight="bold", fontsize=10.3, pad=6)
    ax_d.grid(axis="x", linestyle=":", alpha=0.35)

    save_v2(fig, "Extended_Data_Fig4_Drug_Resistance_Architecture_V2")


def main():
    make_ext_fig1_v2()
    make_ext_fig2_v2()
    make_ext_fig3_v2()
    make_ext_fig4_v2()


if __name__ == "__main__":
    main()
