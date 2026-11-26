# -*- coding: utf-8 -*-
"""
Fig. 6 (Nature Editorial V2): Convergent positive selection & low-cost DR dual engine.
  a. Circos-style Circular 4.41 Mb H37Rv Chromosome Track showing all 229 candidates,
     18 filtered clonal hitchhiking artifacts, and 9 Tier-1 cross-lineage convergent genes
  b. Multi-dimensional Bubble-Matrix (color = max Fst, size = Fst/support, functional-axis ribbon)
  c. Paired Prevalence Dumbbell (Rooted -> Global) + Log-Scale Odds Ratio Forest Plot (6 low-cost DR mutations)
  d. Editorial Vector Infographic of the Two-Phase Evolutionary Model
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch, Circle, Wedge, Patch

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, PAL_GROUP, PAL_AXIS, CMAP_BORDEAUX,
    set_nature_v2_style, save_v2,
)


def main():
    print("[Fig. 6 V2] Generating Circos Chromosome + Bubble-Matrix + Dumbbell/Forest + Editorial Schematic...")
    set_nature_v2_style()

    df_coord = pd.read_csv(os.path.join(BASE_DIR, "Fig. 6", "candidate_gene_coordinates.tsv"), sep="\t")
    df_int = pd.read_csv(os.path.join(BASE_DIR, "Fig. 6", "candidate_gene_integrated_Fst_pi_TajimaD_evidence.tsv"), sep="\t")
    df_pur = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S3_Purified_Convergent_Adaptive_Genes.tsv"), sep="\t")
    df_art = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S3b_Flagged_Clade_Hitchhiking_Candidates.tsv"), sep="\t")
    df_drg = pd.read_csv(os.path.join(BASE_DIR, "Extended Data Fig. 4", "DR_mutation_global_enriched.tsv"), sep="\t")

    gene_max_fst = df_int.groupby("gene_name")["max_Fst"].max().reset_index()
    df_man = df_coord[["gene_name", "gene_mid", "n_mutations", "n_comparisons"]].drop_duplicates("gene_name").merge(
        gene_max_fst, on="gene_name", how="left"
    )
    t1_set = set(df_pur[df_pur["Purified_Tier"] == "Tier1_Gold_CrossLineage_MultiHit"]["gene_name"])
    t2_set = set(df_pur[df_pur["Purified_Tier"] == "Tier2_WithinLineage_MultiHit"]["gene_name"])
    art_set = set(df_art["gene_name"])

    axis_tags = {
        "aftD": "Envelope", "pks12": "Envelope", "pks6": "Envelope",
        "ceoB": "Cation", "Rv3236c": "Cation",
        "pknA": "Kinase", "pknJ": "Kinase",
        "Rv0095c": "Adaptive", "Rv0104": "Adaptive",
    }

    fig = plt.figure(figsize=(16.5, 12.2))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.15, 1.0], width_ratios=[1.05, 1.15], hspace=0.28, wspace=0.24)

    # ==========================================================================
    # PANEL a: Circos-Style Circular 4.41 Mb MTBC Chromosome Selection Track
    # ==========================================================================
    ax_a = fig.add_subplot(gs[0, 0], projection="polar")
    genome_len = 4411532.0
    ax_a.set_theta_zero_location("N")
    ax_a.set_theta_direction(-1)

    # Track bands: inner n_comp track [0.42..0.60], main Fst track [0.65..1.02], outer genome ring [1.06..1.10]
    th_full = np.linspace(0, 2 * np.pi, 360)
    ax_a.fill_between(th_full, 0.65, 1.02, color="#FAF8F4", zorder=0)
    ax_a.fill_between(th_full, 0.42, 0.60, color="#F3F6F8", zorder=0)
    for r_g in [0.65, 0.835, 1.02]:
        ax_a.plot(th_full, np.full_like(th_full, r_g), color="#D6D3CD", lw=0.55, ls=":" if r_g == 0.835 else "-", zorder=1)
    # Outer chromosome ideogram ring
    ax_a.plot(th_full, np.full_like(th_full, 1.07), color="#2B2D42", lw=3.2, zorder=2)

    # Mb ticks around chromosome
    for mb in [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]:
        ang = 2 * np.pi * (mb * 1e6 / genome_len)
        ax_a.plot([ang, ang], [1.07, 1.115], color="#2B2D42", lw=0.8)
        ax_a.text(ang, 1.165, f"{mb:.1f} Mb" if mb == int(mb) else f"{mb:.1f}", ha="center", va="center",
                  fontsize=6.6, color="#525B66")

    # Inner ring bars = n_comparisons
    for _, r in df_man.iterrows():
        ang = 2 * np.pi * (r["gene_mid"] / genome_len)
        gn = r["gene_name"]
        h_in = 0.42 + (min(r["n_comparisons"], 5) / 5.0) * 0.17
        if gn in t1_set:
            c = PAL_AXIS[axis_tags[gn]]
            ax_a.plot([ang, ang], [0.42, h_in], color=c, lw=1.6, zorder=4)
        elif gn in t2_set:
            ax_a.plot([ang, ang], [0.42, h_in], color="#2A6F82", lw=1.0, alpha=0.8, zorder=3)
        else:
            ax_a.plot([ang, ang], [0.42, h_in], color="#CBD5E0", lw=0.45, alpha=0.6, zorder=2)

    # Middle ring scatter = max_Fst
    bg = df_man[~df_man["gene_name"].isin(t1_set | t2_set | art_set)]
    art = df_man[df_man["gene_name"].isin(art_set)]
    t2 = df_man[df_man["gene_name"].isin(t2_set)]
    t1 = df_man[df_man["gene_name"].isin(t1_set)]

    def to_r(fst_arr):
        return 0.65 + np.clip(fst_arr, 0, 1.0) * 0.36

    ax_a.scatter(2 * np.pi * bg["gene_mid"] / genome_len, to_r(bg["max_Fst"]),
                 s=11, color="#BFC6CE", alpha=0.65, linewidths=0, zorder=2)
    ax_a.scatter(2 * np.pi * art["gene_mid"] / genome_len, to_r(art["max_Fst"]),
                 s=24, marker="x", color="#6B7280", linewidths=0.85, zorder=3)
    ax_a.scatter(2 * np.pi * t2["gene_mid"] / genome_len, to_r(t2["max_Fst"]),
                 s=34, color="#2A6F82", edgecolor="white", lw=0.5, zorder=4)

    # Tier 1 genes with clean non-overlapping outer radial callouts
    t1_sorted = t1.sort_values("gene_mid").reset_index(drop=True)
    # Manually nudge clustered genes (early Mb & ~2.3 Mb pks12/pknJ) cleanly around the arc
    ang_override_deg = {
        "pknA": 4.0, "Rv0095c": 16.0, "Rv0104": 28.0, "aftD": 40.0, "pks6": 53.0,
        "pks12": 182.0, "pknJ": 200.0,
    }
    for _, r in t1_sorted.iterrows():
        gn = r["gene_name"]
        ang = 2 * np.pi * (r["gene_mid"] / genome_len)
        rv = to_r(r["max_Fst"])
        c = PAL_AXIS[axis_tags[gn]]
        ax_a.scatter([ang], [rv], s=68, color=c, edgecolor="white", lw=0.8, zorder=6)
        ang_lbl = np.deg2rad(ang_override_deg[gn]) if gn in ang_override_deg else ang
        ax_a.plot([ang, ang_lbl], [rv, 1.25], color=c, lw=0.85, alpha=0.85, zorder=5)
        deg = np.rad2deg(ang_lbl)
        ha = "left" if (0 <= deg <= 180) else "right"
        ax_a.text(ang_lbl, 1.29, f"{gn} ({int(r['n_mutations'])}m)", ha=ha, va="center",
                  fontsize=7.8, fontweight="bold", fontstyle="italic", color=c)

    ax_a.text(0, 0, "H37Rv\n4.41 Mb\n\n9 Tier-1\nCross-L2+L4\nGenes\n(18 artifacts\nfiltered)",
              ha="center", va="center", fontsize=7.6, fontweight="bold", color="#2B2D42")
    ax_a.set_ylim(0, 1.50)
    ax_a.axis("off")
    ax_a.set_title("a   Circular MTBC genome scan & clonal hitchhiking filtration",
                   loc="left", fontweight="bold", fontsize=10.3, pad=6)

    # ==========================================================================
    # PANEL b: Multi-Dimensional Editorial Bubble-Matrix (Tier 1 & Tier 2)
    # ==========================================================================
    ax_b = fig.add_subplot(gs[0, 1])
    t12 = df_pur[df_pur["Purified_Tier"].isin(["Tier1_Gold_CrossLineage_MultiHit", "Tier2_WithinLineage_MultiHit"])].copy()
    axis_order = ["aftD", "pks12", "pks6", "ceoB", "Rv3236c", "pknA", "pknJ", "Rv0095c", "Rv0104",
                  "Rv3521", "Rv3190c", "argS", "pheT", "Rv2800"]
    t12["ord"] = t12["gene_name"].apply(lambda g: axis_order.index(g) if g in axis_order else 99)
    t12 = t12.sort_values("ord").reset_index(drop=True)

    fst_cols = ["Fst_global_vs_rooted", "Fst_L22_vs_L21", "Fst_L41", "Fst_L43", "Fst_L44", "Fst_L48"]
    col_lbls = ["Global vs\nRooted", "L2.2 vs\nL2.1", "L4.1 vs\nother L4", "L4.3 vs\nother L4", "L4.4 vs\nother L4", "L4.8 vs\nother L4"]

    n_rows, n_cols = len(t12), len(fst_cols)
    # Light alternating row backgrounds
    for i in range(n_rows):
        if i % 2 == 1:
            ax_b.axhspan(i - 0.5, i + 0.5, color="#FAF8F4", zorder=0)
    ax_b.axhline(8.5, color="#2B2D42", ls="--", lw=1.2, zorder=2)

    # Left colored physiological-axis strip
    for i, (_, r) in enumerate(t12.iterrows()):
        tag = axis_tags.get(r["gene_name"], "Tier 2")
        ax_b.scatter([-0.68], [i], s=95, marker="s", color=PAL_AXIS[tag], edgecolor="none", zorder=3)
        for j, fc in enumerate(fst_cols):
            val = float(r[fc]) if pd.notna(r[fc]) else 0.0
            sz = 28 + 260 * (val ** 1.35)
            col_v = CMAP_BORDEAUX(val)
            ax_b.scatter([j], [i], s=sz, color=col_v, edgecolor="#333840" if val >= 0.55 else "#CBD5E0",
                         linewidth=0.65 if val >= 0.55 else 0.35, zorder=4)
            if val >= 0.55:
                ax_b.text(j, i, f"{val:.2f}", ha="center", va="center",
                          fontsize=6.5, fontweight="bold", color="white" if val >= 0.70 else "#1E2229", zorder=5)

    ax_b.set_xlim(-1.05, n_cols - 0.35)
    ax_b.set_ylim(n_rows - 0.5, -0.5)
    ax_b.set_xticks(range(n_cols))
    ax_b.set_xticklabels(col_lbls, fontsize=7.8)
    ax_b.set_yticks(range(n_rows))
    ax_b.set_yticklabels([f"{r['gene_name']}  ({int(r['n_mutations'])} mut)" for _, r in t12.iterrows()],
                         fontsize=8.0, fontstyle="italic")
    ax_b.spines["left"].set_visible(False)
    ax_b.tick_params(axis="y", length=0)

    # Colorbar for Fst
    sm = plt.cm.ScalarMappable(cmap=CMAP_BORDEAUX, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=ax_b, shrink=0.62, pad=0.02, aspect=22)
    cb.set_label("Gene max $F_{ST}$ (size & color)", fontsize=8.0)
    cb.ax.tick_params(labelsize=7.5)

    # Axis legend above/below
    h_ax = [Patch(facecolor=PAL_AXIS[k], label=k) for k in ["Envelope", "Cation", "Kinase", "Adaptive", "Tier 2"]]
    ax_b.legend(handles=h_ax, loc="upper center", bbox_to_anchor=(0.48, -0.11), ncol=5,
                frameon=False, fontsize=7.5, handlelength=1.0, columnspacing=0.9)
    ax_b.set_title("b   Purified cross-lineage (Tier 1) & within-L4 (Tier 2) convergence matrix",
                   loc="left", fontweight="bold", fontsize=10.3, pad=6)

    # ==========================================================================
    # PANEL c: Split Dumbbell (Prevalence) + Gradient Forest Plot (OR & 95% CI)
    # ==========================================================================
    gs_c = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 0], width_ratios=[1.0, 1.35], wspace=0.08)
    ax_c1 = fig.add_subplot(gs_c[0])
    ax_c2 = fig.add_subplot(gs_c[1])

    df_drg = df_drg.sort_values("odds_ratio", ascending=True).reset_index(drop=True)
    y_c = np.arange(len(df_drg))
    n_g, n_r = 15257, 941

    # Left subpanel: Rooted -> Global paired frequency dumbbell
    for yi, (_, r) in zip(y_c, df_drg.iterrows()):
        fg, fr = r["freq_global"] * 100, r["freq_rooted"] * 100
        ax_c1.plot([fr, fg], [yi, yi], color="#C44536", lw=2.2, alpha=0.65, zorder=2)
        ax_c1.scatter([fr], [yi], s=48, color=PAL_GROUP["Rooted"], edgecolor="white", lw=0.7, zorder=4)
        ax_c1.scatter([fg], [yi], s=58, color=PAL_GROUP["Global"], edgecolor="white", lw=0.7, zorder=4)
    drug_lbls = []
    for _, r in df_drg.iterrows():
        ds = str(r["drugs"]).replace("ethionamide;isoniazid", "INH/ETH").replace("streptomycin", "STR").replace("ethambutol", "EMB").replace("isoniazid", "INH").replace("rifampicin", "RIF")
        drug_lbls.append(f"{r['gene_name']} {r['change']}\n({ds})")
    ax_c1.set_yticks(y_c)
    ax_c1.set_yticklabels(drug_lbls, fontsize=7.7)
    ax_c1.set_xlim(0, 33)
    ax_c1.set_xlabel("Prevalence (%)")
    ax_c1.grid(axis="x", linestyle=":", alpha=0.35)
    h_db = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=PAL_GROUP["Rooted"], markersize=6.5, label="Rooted"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=PAL_GROUP["Global"], markersize=6.5, label="Global"),
    ]
    ax_c1.legend(handles=h_db, loc="lower right", frameon=True, facecolor="#FAF8F4", edgecolor="#D6D3CD", fontsize=7.0)
    ax_c1.set_title("c   Low-fitness-cost DR enrichment", loc="left", fontweight="bold", fontsize=10.3, pad=6)

    # Right subpanel: Log-scale OR Forest Plot
    ax_c2.axvline(1.0, color="#6B7280", ls="--", lw=0.9, zorder=1)
    for yi, (_, r) in zip(y_c, df_drg.iterrows()):
        a, c = float(r["global"]), float(r["rooted"])
        or_v = float(r["odds_ratio"])
        se = np.sqrt(1/a + 1/(n_g - a) + 1/c + 1/(n_r - c))
        lo, hi = np.exp(np.log(or_v) - 1.96 * se), np.exp(np.log(or_v) + 1.96 * se)
        ax_c2.plot([lo, hi], [yi, yi], color="#2B2D42", lw=1.6, solid_capstyle="round", zorder=2)
        ax_c2.scatter([or_v], [yi], s=70, color="#C44536", edgecolor="white", lw=0.8, zorder=4)
        fdr_s = f"{r['FDR']:.1e}" if r["FDR"] < 1e-3 else f"{r['FDR']:.3f}"
        ax_c2.text(hi * 1.09, yi, f"OR {or_v:.2f}\nFDR {fdr_s}", va="center", fontsize=7.0, color="#2B2D42")
    ax_c2.set_xscale("log")
    ax_c2.set_xlim(0.88, 16.5)
    ax_c2.set_xticks([1, 2, 4, 8])
    ax_c2.set_xticklabels(["1", "2", "4", "8"])
    ax_c2.set_yticks(y_c)
    ax_c2.set_yticklabels([])
    ax_c2.spines["left"].set_visible(False)
    ax_c2.tick_params(axis="y", length=0)
    ax_c2.set_xlabel("Odds ratio (Global vs. Rooted, 95% CI)")
    ax_c2.grid(axis="x", linestyle=":", alpha=0.35)

    # ==========================================================================
    # PANEL d: Editorial Vector Infographic of Two-Phase Evolutionary Model
    # ==========================================================================
    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.set_xlim(0, 10)
    ax_d.set_ylim(0, 10)
    ax_d.axis("off")
    ax_d.set_title("d   Two-phase evolutionary architecture of cosmopolitan TB", loc="left", fontweight="bold", fontsize=10.3, pad=6)

    def card(x, y, w, h, title, body, fc, ec, bar_c, tcol):
        ax_d.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.16,rounding_size=0.25",
                                      facecolor=fc, edgecolor=ec, lw=1.2, zorder=2))
        # Left editorial color accent bar
        ax_d.add_patch(FancyBboxPatch((x - 0.08, y + 0.08), 0.18, h - 0.16, boxstyle="round,pad=0.02,rounding_size=0.08",
                                      facecolor=bar_c, edgecolor="none", zorder=3))
        ax_d.text(x + 0.30, y + h - 0.25, title, ha="left", va="top", fontsize=8.5, fontweight="bold", color=tcol, zorder=4)
        ax_d.text(x + 0.30, y + h/2 - 0.24, body, ha="left", va="center", fontsize=7.5, color="#2B2D42", linespacing=1.35, zorder=4)

    card(0.25, 6.95, 9.5, 2.65,
         "PHASE 1 · Pre-Antibiotic Convergent Adaptation (Independent in L2 & L4)",
         "• Cell-envelope remodeling:  aftD (arabinan/LAM), pks12 (MPM/CD1c), pks6\n"
         "• Intraphagosomal cation homeostasis:  ceoB (Trk K+ uptake), Rv3236c\n"
         "• Upstream STPK signaling:  pknA, pknJ   +   conserved loci (Rv0095c, Rv0104)",
         fc="#FDF7F5", ec="#E4B8AE", bar_c="#C44536", tcol="#8C261D")

    ax_d.annotate("", xy=(5.0, 5.75), xytext=(5.0, 6.72), arrowprops=dict(arrowstyle="-|>", color="#2B2D42", lw=1.8))
    ax_d.text(5.25, 6.24, "Broad host-range competence (85 vs 46 countries, P = 0.0033)",
              va="center", fontsize=7.4, fontweight="bold", color="#2B2D42")

    card(0.25, 3.05, 9.5, 2.48,
         "PHASE 2 · Mid-20th-Century Chemotherapy Bottleneck (Post-1950)",
         "• Global clades enrich low-fitness-cost, transmission-permissive alleles:\n"
         "   rpsL K43R (OR 3.57; 10.95× in L2.2 vs L2.1), embB M306V, katG S315T, rpoB S450L\n"
         "• Total MDR rate alone fails to explain success (rooted L2.1: 44.3% MDR yet collapsed)",
         fc="#F3F8FA", ec="#B4D2DC", bar_c="#2A6F82", tcol="#1D4F5E")

    ax_d.annotate("", xy=(2.5, 2.15), xytext=(4.2, 2.85), arrowprops=dict(arrowstyle="-|>", color="#C44536", lw=1.7))
    ax_d.annotate("", xy=(7.5, 2.15), xytext=(5.8, 2.85), arrowprops=dict(arrowstyle="-|>", color="#6B7280", lw=1.7))

    card(0.25, 0.20, 4.5, 1.78,
         "Cosmopolitan Generalists",
         "L2.2, L4.1, L4.3, L4.4, L4.8\n98.0% 21st-century Ne retention\n(modern demographic resilience)",
         fc="#FBECE8", ec="#C44536", bar_c="#C44536", tcol="#8C261D")

    card(5.25, 0.20, 4.5, 1.78,
         "Regionally Rooted Specialists",
         "L2.1, L4.5 (and localized L3 clades)\n>91% peak collapse (29.1% since 2000)\n(contraction under chemotherapy)",
         fc="#F2F4F6", ec="#8D99AE", bar_c="#6B7280", tcol="#374151")

    save_v2(fig, "Fig6_Convergent_Adaptation_and_LowCost_DR_V2")


if __name__ == "__main__":
    main()
