# -*- coding: utf-8 -*-
"""
Fig. 5 (Nature Editorial V2): Sublineage-resolved Bayesian phylogeography.
  a-f. Six flagship Circular / Fan Time-Calibrated BEAST MCC Trees with concentric century
       rings and outer geographic tip bands contrasting rooted vs cosmopolitan radiations
  g.   Time-calibrated tMRCA & 95% HPD Caterpillar Plot across all 14 sublineages
  h.   Ancestral Origin Posterior Probability Stacked Bar Synthesis across all 14 sublineages
"""

import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Patch

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, CLEAN_DIR, PAL_GROUP, PAL_REGION,
    set_nature_v2_style, save_v2, draw_circular_beast_tree,
)


def parse_root_full_probs(tree_path):
    txt = open(tree_path, encoding="utf-8", errors="ignore").read()
    m = re.search(r"tree\s+\S+\s*=\s*(?:\[&R\]\s*)?(.+);", txt, re.S)
    annots = re.findall(r"\[&([^\]]+)\]", m.group(1))
    root_a = annots[-1]

    def gval(s, k):
        mm = re.search(r"(?:^|,)" + re.escape(k) + r"=(\{[^}]+\}|\"[^\"]+\"|[^,]+)", s)
        return mm.group(1) if mm else ""

    lset = gval(root_a, "location.set").strip("{}").replace('"', "").split(",")
    lprob = [float(x) for x in gval(root_a, "location.set.prob").strip("{}").split(",")]
    return dict(zip(lset, lprob))


def main():
    print("[Fig. 5 V2] Generating 6 Circular Fan MCC Trees + tMRCA Caterpillar + Origin Posterior Matrix...")
    set_nature_v2_style()

    df_qc = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S4_BEAST_Tree_QC_and_Outlier_Audit.tsv"), sep="\t").set_index("sublineage")
    df_s2 = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv"), sep="\t").set_index("sublineage")

    # 6 flagship circular trees showcasing sharpest ecological contrasts & high-confidence origins
    flagship = [
        ("lineage2.1", "a  L2.1 (Rooted)"),
        ("lineage2.2", "b  L2.2 (Global)"),
        ("lineage4.5", "c  L4.5 (Rooted)"),
        ("lineage4.8", "d  L4.8 (Global)"),
        ("lineage4.7", "e  L4.7 (Europe)"),
        ("lineage1.1", "f  L1.1 (Ancient L1)"),
    ]

    fig = plt.figure(figsize=(16.5, 11.8))
    outer = gridspec.GridSpec(2, 1, height_ratios=[1.22, 1.0], hspace=0.28)

    # Top row: 2x3 or 1x6? 2x3 gives HUGE gorgeous crisp circular trees! Wait: 2 rows x 3 cols in top half, OR 1 row x 6?
    # 2x3 circular trees in top 60% of figure is MUCH bigger and cleaner than 1x6!
    gs_top = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=outer[0], hspace=0.22, wspace=0.12)

    for idx, (sn, title_lbl) in enumerate(flagship):
        ax = fig.add_subplot(gs_top[idx // 3, idx % 3])
        tpath = os.path.join(BASE_DIR, "Fig. 5", f"{sn}_location_tree_with_trait.tree")
        qc = df_qc.loc[sn]
        grp = df_s2.loc[sn, "ecological_group"]
        loc_pretty = str(qc["root_top1_loc"]).replace("_", " ")
        badge = f"Origin: {loc_pretty} (PP = {qc['root_top1_PP']:.2f})  |  tMRCA = {qc['root_height_yr']:.0f} y"
        draw_circular_beast_tree(ax, tpath, title_lbl, badge, PAL_GROUP[grp])

    # Bottom row: g (tMRCA Caterpillar Plot) + h (Ancestral Origin Posterior Stacked Bars)
    gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], width_ratios=[1.0, 1.32], wspace=0.26)

    # Order all 14 sublineages by ecological group & tMRCA
    order_14 = [
        "lineage1.1", "lineage1.2", "lineage3.1",
        "lineage2.2", "lineage4.1", "lineage4.4", "lineage4.3", "lineage4.8",
        "lineage2.1", "lineage4.5",
        "lineage4.6", "lineage4.2", "lineage4.7", "lineage4.9",
    ][::-1]

    # --- Panel g: Time-Calibrated tMRCA & 95% HPD Caterpillar Plot ---
    ax_g = fig.add_subplot(gs_bot[0])
    y_pos = np.arange(len(order_14))
    for yi, sn in enumerate(order_14):
        qc = df_qc.loc[sn]
        grp = df_s2.loc[sn, "ecological_group"]
        col = PAL_GROUP[grp]
        h_val = float(qc["root_height_yr"])
        lo, hi = [float(x) for x in str(qc["root_height_95HPD"]).split("-")]
        ax_g.plot([lo, hi], [yi, yi], color=col, lw=2.3, alpha=0.65, solid_capstyle="round", zorder=2)
        ax_g.scatter([h_val], [yi], s=58, color=col, edgecolor="white", lw=0.8, zorder=4)
        ax_g.text(hi + 26, yi, f"{h_val:.0f} y", va="center", fontsize=7.0, color="#2B2D42")

    ax_g.set_yticks(y_pos)
    ax_g.set_yticklabels([f"{s.replace('lineage','L')} [{df_s2.loc[s,'ecological_group'][0]}]" for s in order_14], fontsize=7.9, fontweight="bold")
    ax_g.set_xlim(250, 1950)
    ax_g.set_xlabel("Estimated tMRCA (years before 2025; 95% HPD)")
    ax_g.set_title("g   Divergence depths (tMRCA & 95% HPD) across 14 sublineages", loc="left", fontweight="bold", fontsize=10.2, pad=6)
    ax_g.grid(axis="x", linestyle=":", alpha=0.35)

    # --- Panel h: Ancestral Origin Posterior Probability Stacked Bar Chart ---
    ax_h = fig.add_subplot(gs_bot[1])
    reg_order = [
        "East_Asia", "Southeast_Asia", "South_Asia", "Western_Asia", "Central_Asia",
        "Europe", "Africa", "Northern_America", "South_America", "Oceania",
    ]
    for yi, sn in enumerate(order_14):
        if sn == "lineage4.1":
            tpath = os.path.join(CLEAN_DIR, "lineage4.1_location_tree_with_trait_clean499.tree")
        else:
            tpath = os.path.join(BASE_DIR, "Fig. 5", f"{sn}_location_tree_with_trait.tree")
        probs = parse_root_full_probs(tpath)
        left = 0.0
        for reg in reg_order:
            p = probs.get(reg, 0.0) * 100.0
            if p > 0.3:
                ax_h.barh(yi, p, left=left, height=0.66, color=PAL_REGION[reg], edgecolor="white", lw=0.35)
            left += p
        top_pp = float(df_qc.loc[sn, "root_top1_PP"])
        top_loc = str(df_qc.loc[sn, "root_top1_loc"]).replace("_", " ")
        tag = f"● {top_loc} ({top_pp:.2f})" if top_pp >= 0.65 else f"Uncertain ({top_loc} {top_pp:.2f})"
        tcol = "#1E2229" if top_pp >= 0.65 else "#6B7280"
        ax_h.text(101.8, yi, tag, va="center", fontsize=7.2, fontweight="bold" if top_pp >= 0.65 else "normal", color=tcol)

    ax_h.axvline(65, color="#2B2D42", ls="--", lw=0.8, alpha=0.6)
    ax_h.set_yticks(y_pos)
    ax_h.set_yticklabels([s.replace("lineage", "L") for s in order_14], fontsize=7.9)
    ax_h.set_xlim(0, 152)
    ax_h.set_xticks([0, 25, 50, 65, 75, 100])
    ax_h.set_xlabel("Ancestral root location posterior probability (%; dashed line: PP = 0.65 threshold)")
    ax_h.set_title("h   Ancestral origin posterior distributions (● high-confidence PP ≥ 0.65)", loc="left", fontweight="bold", fontsize=10.2, pad=6)

    handles = [Patch(facecolor=PAL_REGION[r], label=r.replace("_", " ")) for r in reg_order]
    ax_h.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=5,
                frameon=True, facecolor="#FAF8F4", edgecolor="#D6D3CD", fontsize=7.5, columnspacing=0.9)

    save_v2(fig, "Fig5_Sublineage_Phylogeography_Trees_V2")


if __name__ == "__main__":
    main()
