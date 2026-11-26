# -*- coding: utf-8 -*-
"""
Fig. 1 (Nature Editorial V2): Century-spanning global genomic atlas of MTBC.
  a. Robinson-projected global map of 107,309 geolocated isolates (120 countries)
     with editorial donut-pie markers & key-country callouts
  b. Concentric Dual-Ring Sunburst (inner: taxonomic lineages; outer: WHO DR spectrum; n = 162,869)
  c. Smooth Stacked Area Trajectory (1995-2025) with public-health milestones
  d. Marimekko Proportional-Width Mosaic Plot of regional lineage architecture
"""

import os
import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import Wedge, Circle, Patch, FancyBboxPatch

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, PAL_LINEAGE,
    set_nature_v2_style, save_v2, robinson_project, draw_robinson_basemap,
)

MAJOR_LIN = ["lineage4", "lineage2", "lineage1", "lineage3", "La1", "lineage6", "lineage5", "Others"]
MAJOR_LBL = {
    "lineage4": "L4 (Euro-American)",
    "lineage2": "L2 (East Asian)",
    "lineage1": "L1 (Indo-Oceanic)",
    "lineage3": "L3 (East African-Indian)",
    "La1": "La1 (M. bovis)",
    "lineage6": "L6 (West African 2)",
    "lineage5": "L5 (West African 1)",
    "Others": "Others (L7-9, La2-3, M. can.)",
}


def draw_donut_pie(ax, x, y, fracs, colors, r_outer, hole_ratio=0.38):
    ang = 90.0
    for frac, col in zip(fracs, colors):
        if frac <= 0.002:
            continue
        ang2 = ang - frac * 360.0
        w = Wedge((x, y), r_outer, ang2, ang, width=r_outer * (1.0 - hole_ratio),
                  facecolor=col, edgecolor="white", linewidth=0.35, zorder=5)
        ax.add_patch(w)
        ang = ang2
    inner = Circle((x, y), r_outer * hole_ratio, facecolor="#FFFFFF", edgecolor="none", zorder=5.2, alpha=0.92)
    ax.add_patch(inner)


def main():
    print("[Fig. 1 V2] Generating Nature-editorial Fig1_Global_MTBC_Atlas...")
    set_nature_v2_style()

    df_geo = pd.read_csv(os.path.join(BASE_DIR, "Fig. 1", "data.csv")).fillna(0)
    df_s1b = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S1b_MainLineage_and_DR_Breakdown.tsv"), sep="\t")
    df_lin = df_s1b[df_s1b["Category"] == "Main_Lineage"].copy()
    df_dr = df_s1b[df_s1b["Category"] == "Drug_Resistance_Type"].copy()

    wb = openpyxl.load_workbook(os.path.join(BASE_DIR, "Fig. 1", "data.xlsx"), read_only=True, data_only=True)
    r8 = list(wb["Sheet8"].iter_rows(values_only=True))
    df_yr = pd.DataFrame(r8[1:], columns=r8[0])
    r15 = list(wb["Sheet15"].iter_rows(values_only=True))
    df_reg = pd.DataFrame(r15[1:], columns=r15[0])
    wb.close()

    lin_all = ["La1", "La2", "La3", "lineage1", "lineage2", "lineage3", "lineage4",
               "lineage5", "lineage6", "lineage7", "lineage8", "lineage9", "M.canetti"]
    df_yr["Year"] = pd.to_numeric(df_yr["Year"], errors="coerce")
    df_yr = df_yr.dropna(subset=["Year"]).copy()
    df_yr["Year"] = df_yr["Year"].astype(int)
    for c in lin_all:
        df_yr[c] = pd.to_numeric(df_yr[c], errors="coerce").fillna(0)
        df_reg[c] = pd.to_numeric(df_reg[c], errors="coerce").fillna(0)
    df_yr["Others"] = df_yr[["La2", "La3", "lineage7", "lineage8", "lineage9", "M.canetti"]].sum(axis=1)
    df_reg["Others"] = df_reg[["La2", "La3", "lineage7", "lineage8", "lineage9", "M.canetti"]].sum(axis=1)

    fig = plt.figure(figsize=(16.0, 10.2))
    gs = gridspec.GridSpec(2, 3, height_ratios=[1.48, 1.0], width_ratios=[1.05, 1.12, 1.15], hspace=0.24, wspace=0.28)

    # ==========================================================================
    # PANEL a: Robinson-Projected World Map with Editorial Donut Pies
    # ==========================================================================
    ax_a = fig.add_subplot(gs[0, :])
    draw_robinson_basemap(ax_a, land_fc="#ECEAE2", land_ec="#A39E93", land_lw=0.35, ocean_fc="#F3F7F9")

    other_cols = ["La2", "La3", "lineage7", "lineage8", "lineage9", "M.canetti"]
    df_geo["Others"] = df_geo[other_cols].sum(axis=1)
    df_geo_s = df_geo.sort_values("Size", ascending=True)
    pie_cols = [PAL_LINEAGE[c] for c in MAJOR_LIN]

    for _, row in df_geo_s.iterrows():
        sz = float(row["Size"])
        lon, lat = float(row["Longitude"]), float(row["Latitude"])
        if sz <= 0 or lat < -56:
            continue
        vals = np.array([float(row[c]) for c in MAJOR_LIN])
        tot = vals.sum()
        if tot <= 0:
            continue
        fracs = vals / tot
        rx, ry = robinson_project(lon, lat)
        # Robinson x spans ~[-2.66, 2.66], y spans ~[-0.95, 1.32]
        rad = 0.025 + 0.026 * np.log10(max(sz, 1.0))
        if sz >= 45:
            draw_donut_pie(ax_a, rx, ry, fracs, pie_cols, rad, hole_ratio=0.36)
        else:
            dom = pie_cols[int(np.argmax(fracs))]
            ax_a.scatter(rx, ry, s=16, color=dom, edgecolor="white", linewidth=0.35, zorder=4)

    # Subtle key-country micro labels (top hubs)
    hub_labels = [
        ("United States", -0.35, 0.14, "USA (11.5k)"),
        ("China", 0.34, 0.13, "China (11.0k)"),
        ("South Africa", 0.30, -0.10, "S. Africa (9.2k)"),
        ("India", 0.25, -0.14, "India (7.1k)"),
        ("Australia", 0.30, -0.14, "Australia (6.1k)"),
        ("Peru", -0.34, -0.08, "Peru (4.4k)"),
    ]
    geo_idx = df_geo.set_index("Country")
    for ctry, dx, dy, lbl in hub_labels:
        if ctry in geo_idx.index:
            rx, ry = robinson_project(geo_idx.loc[ctry, "Longitude"], geo_idx.loc[ctry, "Latitude"])
            ax_a.annotate(lbl, xy=(rx, ry), xytext=(rx + dx, ry + dy),
                          fontsize=7.3, fontweight="bold", color="#2B2D42", zorder=7,
                          bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="#CBD5E0", alpha=0.88),
                          arrowprops=dict(arrowstyle="-", color="#64748B", lw=0.65))

    # Sleek size scale in bottom-left Pacific
    sx0, sy0 = -2.35, -0.58
    ax_a.text(sx0, sy0 + 0.28, "Country sample size", fontsize=8.0, fontweight="bold", color="#2B2D42")
    for n_ex, y_off in [(100, 0.16), (1000, 0.02), (10000, -0.15)]:
        r_ex = 0.025 + 0.026 * np.log10(n_ex)
        ax_a.add_patch(Circle((sx0 + 0.10, sy0 + y_off), r_ex, facecolor="none", edgecolor="#4A5568", lw=0.8, zorder=6))
        ax_a.text(sx0 + 0.27, sy0 + y_off, f"n = {n_ex:,}", va="center", fontsize=7.5, color="#333840")

    handles = [Patch(facecolor=PAL_LINEAGE[c], edgecolor="white", label=MAJOR_LBL[c]) for c in MAJOR_LIN]
    ax_a.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.54, -0.05), ncol=4,
                frameon=True, facecolor="#FAF8F4", edgecolor="#D6D3CD", fontsize=7.8, columnspacing=1.0)
    ax_a.set_title("a   Global distribution of 107,309 geolocated MTBC isolates across 120 countries (Robinson projection)",
                   loc="left", fontweight="bold", fontsize=10.8, pad=4)

    # ==========================================================================
    # PANEL b: Concentric Dual-Ring Editorial Sunburst (Lineages + WHO DR)
    # ==========================================================================
    ax_b = fig.add_subplot(gs[1, 0])
    ax_b.set_aspect("equal")
    ax_b.axis("off")

    # Inner ring: Lineage proportions across 162,869
    lin_map = dict(zip(df_lin["Label"], df_lin["Count"]))
    other_sum = sum(lin_map.get(k, 0) for k in ["La2", "La3", "lineage7", "lineage9", "M.canetti", "lineage8"])
    inner_vals = [lin_map.get(k, 0) for k in MAJOR_LIN[:-1]] + [other_sum]
    inner_cols = [PAL_LINEAGE[k] for k in MAJOR_LIN]

    # Outer ring: WHO DR spectrum across 162,869
    dr_order = ["Sensitive", "MDR-TB", "Pre-XDR-TB", "XDR-TB", "HR-TB", "RR-TB", "Other"]
    dr_pal = {
        "Sensitive": "#5EA3B8", "MDR-TB": "#C44536", "Pre-XDR-TB": "#7A1C35",
        "XDR-TB": "#380920", "HR-TB": "#D99B26", "RR-TB": "#E59866", "Other": "#C5C2BA",
    }
    dr_map = dict(zip(df_dr["Label"], df_dr["Count"]))
    outer_vals = [dr_map.get(k, 0) for k in dr_order]
    outer_cols = [dr_pal[k] for k in dr_order]

    # Draw inner donut (radius 0.42 -> 0.72) and outer donut (radius 0.76 -> 1.00)
    w_in, _ = ax_b.pie(inner_vals, radius=0.72, startangle=90, counterclock=False,
                       colors=inner_cols, wedgeprops=dict(width=0.28, edgecolor="white", linewidth=0.9))
    w_out, _ = ax_b.pie(outer_vals, radius=1.02, startangle=90, counterclock=False,
                        colors=outer_cols, wedgeprops=dict(width=0.24, edgecolor="white", linewidth=0.9))

    # Center typography
    ax_b.text(0, 0.06, "162,869", ha="center", va="center", fontsize=12.5, fontweight="bold", color="#1E2229")
    ax_b.text(0, -0.10, "MTBC genomes\n1900–2025", ha="center", va="center", fontsize=7.6, color="#525B66")

    # Clean radial labels for major wedges on outer ring
    tot_dr = sum(outer_vals)
    ang = 90.0
    for k, v in zip(dr_order, outer_vals):
        frac = v / tot_dr
        mid = ang - frac * 180.0
        ang -= frac * 360.0
        if frac >= 0.06:
            rad_ang = np.deg2rad(mid)
            lx, ly = 1.18 * np.cos(rad_ang), 1.18 * np.sin(rad_ang)
            ha = "left" if lx >= 0 else "right"
            ax_b.text(lx, ly, f"{k}\n({frac*100:.1f}%)", ha=ha, va="center", fontsize=7.0, fontweight="bold", color="#2B2D42")

    # Inner ring key percentages directly inside large wedges
    tot_in = sum(inner_vals)
    ang = 90.0
    for k, v in zip(MAJOR_LIN, inner_vals):
        frac = v / tot_in
        mid = ang - frac * 180.0
        ang -= frac * 360.0
        if frac >= 0.06:
            rad_ang = np.deg2rad(mid)
            ax_b.text(0.58 * np.cos(rad_ang), 0.58 * np.sin(rad_ang),
                      f"{k.replace('lineage','L')}\n{frac*100:.0f}%",
                      ha="center", va="center", fontsize=6.7, fontweight="bold", color="white")

    ax_b.set_xlim(-1.58, 1.58)
    ax_b.set_ylim(-1.28, 1.28)
    ax_b.set_title("b   Lineage (inner) & WHO DR spectrum (outer)", loc="left", fontweight="bold", fontsize=10.2, pad=4)

    # ==========================================================================
    # PANEL c: Smooth Stacked Area Streamgraph (1995-2025)
    # ==========================================================================
    ax_c = fig.add_subplot(gs[1, 1])
    yr_rec = df_yr[df_yr["Year"] >= 1995].sort_values("Year")
    n_hist = int(df_yr[df_yr["Year"] < 1995][MAJOR_LIN].values.sum())
    ys_stack = [yr_rec[c].values for c in MAJOR_LIN]
    ax_c.stackplot(yr_rec["Year"], ys_stack, colors=[PAL_LINEAGE[c] for c in MAJOR_LIN],
                   alpha=0.88, edgecolor="white", linewidth=0.4)

    # Subtle milestone reference lines
    for yr_m, lbl_m, y_pos in [(2010, "Xpert MTB/RIF\n(2010)", 7800), (2020, "Peak: 2020\n(n=12,736)", 12900)]:
        ax_c.axvline(yr_m, color="#4A5568", linestyle="--", linewidth=0.75, alpha=0.7)
        ax_c.text(yr_m - 0.4, y_pos, lbl_m, ha="right", va="top", fontsize=7.2, color="#2B2D42",
                  bbox=dict(boxstyle="round,pad=0.18", facecolor="#FAF8F4", edgecolor="#D6D3CD", alpha=0.9))

    ax_c.text(0.03, 0.94, f"Historical archive\n(1900–1994): n = {n_hist:,}",
              transform=ax_c.transAxes, va="top", fontsize=7.4, color="#525B66",
              bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="#D6D3CD"))
    ax_c.set_xlim(1995, 2025)
    ax_c.set_ylim(0, 14200)
    ax_c.set_xlabel("Collection year")
    ax_c.set_ylabel("Genomes sequenced per year")
    ax_c.set_title("c   Temporal sequencing trajectory by lineage", loc="left", fontweight="bold", fontsize=10.2, pad=4)
    ax_c.grid(axis="y", linestyle=":", alpha=0.35)

    # ==========================================================================
    # PANEL d: Marimekko Proportional-Width Regional Mosaic Plot
    # ==========================================================================
    ax_d = fig.add_subplot(gs[1, 2])
    reg_agg = df_reg.groupby("Region")[MAJOR_LIN].sum()
    reg_agg = reg_agg[reg_agg.sum(axis=1) > 0]
    reg_tot = reg_agg.sum(axis=1).sort_values(ascending=False)
    reg_agg = reg_agg.loc[reg_tot.index]
    reg_pct = reg_agg.div(reg_tot, axis=0) * 100.0

    grand_tot = reg_tot.sum()
    x_left = 0.0
    gap = 1.4
    x_centers, x_labels = [], []
    for reg in reg_agg.index:
        w = (reg_tot[reg] / grand_tot) * (100.0 - gap * (len(reg_agg) - 1))
        y_bot = 0.0
        for col in MAJOR_LIN:
            h = float(reg_pct.loc[reg, col])
            if h > 0.2:
                ax_d.bar(x_left + w / 2.0, h, width=w, bottom=y_bot,
                         color=PAL_LINEAGE[col], edgecolor="white", linewidth=0.55)
                if h >= 14 and w >= 10:
                    ax_d.text(x_left + w / 2.0, y_bot + h / 2.0,
                              f"{col.replace('lineage','L')}\n{h:.0f}%",
                              ha="center", va="center", fontsize=6.6, fontweight="bold", color="white")
            y_bot += h
        x_centers.append(x_left + w / 2.0)
        x_labels.append(f"{reg}\n({reg_tot[reg]/1000:.1f}k)")
        x_left += w + gap

    ax_d.set_xticks(x_centers)
    ax_d.set_xticklabels(x_labels, fontsize=7.6)
    ax_d.set_xlim(-1, x_left - gap + 1)
    ax_d.set_ylim(0, 103)
    ax_d.set_ylabel("Lineage share (%)")
    ax_d.set_xlabel("Macro-region (bar width scaled to regional sample size)")
    ax_d.set_title("d   Proportional regional mosaic (Marimekko)", loc="left", fontweight="bold", fontsize=10.2, pad=4)

    save_v2(fig, "Fig1_Global_MTBC_Atlas_V2")


if __name__ == "__main__":
    main()
