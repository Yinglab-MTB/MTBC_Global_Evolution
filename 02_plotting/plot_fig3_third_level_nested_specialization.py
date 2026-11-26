# -*- coding: utf-8 -*-
"""
Fig. 3 (Nature Editorial V2): Nested third-level ecological specialization.
  a-f. Six Robinson-projected contrasting sister clades (Terracotta cosmopolitan vs Botanical-Teal localized)
  g.   Circos-style Circular Polar Bar + Shannon Entropy Ring across all 30 third-level sublineages
  h.   Within-clade Nested Divergence Dumbbell Plot (broadest vs narrowest sister in each Global clade)
"""

import os
import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import PolyCollection
from matplotlib.patches import Patch

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, PAL_GROUP, CMAP_TERRACOTTA, CMAP_BOTANICAL,
    set_nature_v2_style, save_v2, robinson_project, draw_robinson_basemap,
    kriging_robinson_surface, load_world_polygons_robinson,
)


def main():
    print("[Fig. 3 V2] Generating Robinson sister contrasts + Circular Polar Ring + Dumbbell plot...")
    set_nature_v2_style()

    df_s3 = pd.read_csv(os.path.join(AUDIT_DIR, "Table_S2c_ThirdLevel_Sublineage_Geography.tsv"), sep="\t")
    s3_idx = df_s3.set_index("sublineage_L3")
    df_c = pd.read_csv(os.path.join(BASE_DIR, "Fig. 1", "data.csv"))
    coord_map = dict(zip(df_c["Country"], zip(df_c["Longitude"], df_c["Latitude"])))

    wb = openpyxl.load_workbook(
        os.path.join(BASE_DIR, "Fig. 3", "二级谱系分布频率.xlsx"), read_only=True, data_only=True
    )

    # 3 flagship cosmopolitan-vs-localized pairs (2x3 Robinson maps)
    featured = [
        ("lineage2.2.1", "a", "Cosmopolitan L2.2 sister", True),
        ("lineage4.1.2", "b", "Cosmopolitan L4.1 sister", True),
        ("lineage4.3.3", "c", "Cosmopolitan L4.3 sister", True),
        ("lineage2.2.2", "d", "Localized L2.2 sister (South Africa)", False),
        ("lineage4.4.2", "e", "Localized L4.4 sister (China)", False),
        ("lineage4.8.1", "f", "Localized L4.8 sister (South Africa)", False),
    ]

    fig = plt.figure(figsize=(16.2, 12.2))
    outer = gridspec.GridSpec(2, 1, height_ratios=[1.12, 1.15], hspace=0.22)
    gs_top = gridspec.GridSpecFromSubplotSpec(2, 3, subplot_spec=outer[0], hspace=0.22, wspace=0.06)

    for i, (sn, letter, tag, is_broad) in enumerate(featured):
        ax = fig.add_subplot(gs_top[i // 3, i % 3])
        ocean = draw_robinson_basemap(ax, land_fc="#F2EFE8", land_ec="#A39E93", land_lw=0.26,
                                      ocean_fc="#F4F8FA", graticule=False, border_ec="#CBD5E0")
        ws = wb[sn]
        rows = list(ws.iter_rows(values_only=True))
        obs_c = {str(r[0]): float(r[1]) for r in rows[1:] if r[0] is not None and isinstance(r[1], (int, float))}
        tot = sum(obs_c.values())

        lons_o, lats_o, vals_o = [], [], []
        for ctry, (lon, lat) in coord_map.items():
            if lat < -56:
                continue
            lons_o.append(lon)
            lats_o.append(lat)
            vals_o.append(obs_c.get(ctry, 0.0) / max(tot, 1.0))

        vmax = max(np.max(vals_o), 0.05)
        RX, RY, grid = kriging_robinson_surface(lons_o, lats_o, vals_o, nx=190, ny=92, range_deg=22.0, nugget=0.06)
        grid_norm = np.clip(grid / vmax, 0.0, 1.0)
        cm = CMAP_TERRACOTTA if is_broad else CMAP_BOTANICAL
        qm = ax.pcolormesh(RX, RY, grid_norm, cmap=cm, vmin=0, vmax=1, shading="gouraud", zorder=2, alpha=0.95)
        qm.set_clip_path(ocean)
        ax.contour(RX, RY, grid_norm, levels=[0.3, 0.6], colors="#4A0E1B" if is_broad else "#0B2F36",
                   linewidths=0.35, alpha=0.45, zorder=2.6)
        ax.add_collection(PolyCollection(load_world_polygons_robinson(), facecolor="none",
                                         edgecolor="#6E6A61", linewidth=0.22, zorder=3))
        nz = np.array(vals_o) > 0
        px, py = robinson_project(np.array(lons_o)[nz], np.array(lats_o)[nz])
        ax.scatter(px, py, s=4.0, color="#1E2229", alpha=0.42, linewidths=0, zorder=3.5)

        meta = s3_idx.loc[sn]
        col = "#C44536" if is_broad else "#1C595A"
        short = sn.replace("lineage", "L")
        ax.text(0.01, 1.02, f" {letter}  {short} ", transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8.2, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.18", facecolor=col, edgecolor="none"))
        ax.text(0.20, 1.02, f"{int(meta['n_countries'])} countries · H={meta['shannon_H']:.2f} | {meta['top1_country']}",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=7.5, color="#2B2D42")

    wb.close()

    # ==========================================================================
    # BOTTOM ROW: g (Circular Polar Ring of 30 L3 sublineages) + h (Nested Dumbbell)
    # ==========================================================================
    gs_bot = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], width_ratios=[1.25, 1.0], wspace=0.24)

    # --- Panel g: Circos-style Circular Polar Bar + Entropy Ring ---
    ax_g = fig.add_subplot(gs_bot[0], projection="polar")
    df_p = df_s3.sort_values(["parent_L2", "n_countries"], ascending=[True, False]).reset_index(drop=True)
    n30 = len(df_p)
    # Leave small visual gaps between parent L2 groups!
    angles, bar_w = [], []
    parents = df_p["parent_L2"].values
    n_groups = len(pd.unique(parents))
    total_slots = n30 + n_groups * 0.8
    slot = 0.0
    group_spans = {}
    for idx, p_name in enumerate(parents):
        if idx > 0 and p_name != parents[idx - 1]:
            slot += 0.8
        ang = 2.0 * np.pi * (slot / total_slots)
        angles.append(ang)
        group_spans.setdefault(p_name, []).append(ang)
        slot += 1.0
    angles = np.array(angles)
    w_rad = 2.0 * np.pi / total_slots * 0.78

    r_base = 52.0
    # Concentric grid circles
    for r_g in [r_base + 25, r_base + 50, r_base + 75]:
        th = np.linspace(0, 2 * np.pi, 200)
        ax_g.plot(th, np.full_like(th, r_g), color="#E2E8F0", lw=0.6, linestyle=":", zorder=0)

    # Bars for country richness
    cols_p = [PAL_GROUP[g] for g in df_p["parent_group"]]
    ax_g.bar(angles, df_p["n_countries"].values, width=w_rad, bottom=r_base,
             color=cols_p, alpha=0.84, edgecolor="white", linewidth=0.5, zorder=3)

    # Inner colored arcs indicating parent 2nd-level clade + crisp non-overlapping labels
    for p_name, ang_list in group_spans.items():
        a0, a1 = min(ang_list) - w_rad * 0.45, max(ang_list) + w_rad * 0.45
        th_arc = np.linspace(a0, a1, 30)
        p_grp = df_p[df_p["parent_L2"] == p_name]["parent_group"].iloc[0]
        ax_g.plot(th_arc, np.full_like(th_arc, r_base - 4.5), color=PAL_GROUP[p_grp], lw=3.4, solid_capstyle="butt", zorder=4)
        amid = 0.5 * (a0 + a1)
        ax_g.text(amid, r_base - 13.5, p_name.replace("lineage", "L"), ha="center", va="center",
                  fontsize=6.5, fontweight="bold", color=PAL_GROUP[p_grp])

    # Outer ring dots = Shannon entropy H (scaled)
    ax_g.scatter(angles, r_base + df_p["shannon_H"].values * 24.0, s=16, color="#1E2229",
                 facecolor="white", linewidth=0.9, zorder=5)

    # True radial tip labels (theta_zero='N', clockwise -> angle deg maps to 90-deg on right, 270-deg on left)
    for ang, (_, row) in zip(angles, df_p.iterrows()):
        deg = np.rad2deg(ang) % 360.0
        if 0 <= deg <= 180:
            rot = 90.0 - deg
            ha = "left"
        else:
            rot = 270.0 - deg
            ha = "right"
        r_txt = r_base + max(row["n_countries"], row["shannon_H"] * 24.0) + 8.5
        ax_g.text(ang, r_txt, row["sublineage_L3"].replace("lineage", "L"),
                  rotation=rot, rotation_mode="anchor", ha=ha, va="center", fontsize=6.5, color="#2B2D42")

    ax_g.set_theta_zero_location("N")
    ax_g.set_theta_direction(-1)
    ax_g.set_ylim(0, r_base + 116)
    ax_g.axis("off")
    ax_g.text(0, 0, "30 L3 Clades\n(bars: countries\ndots: H)", ha="center", va="center",
              fontsize=6.8, fontweight="bold", color="#333840")
    ax_g.set_title("g   Circular breadth & entropy ring across all 30 third-level clades",
                   loc="left", fontweight="bold", fontsize=10.2, pad=8)

    # --- Panel h: Nested Sister-Clade Divergence Dumbbell Plot ---
    ax_h = fig.add_subplot(gs_bot[1])
    # For each parent clade with >=2 L3 sublineages, plot widest vs narrowest sister!
    dumb_rows = []
    for parent, g in df_s3.groupby("parent_L2"):
        if len(g) < 2:
            continue
        g_s = g.sort_values("n_countries", ascending=False)
        hi = g_s.iloc[0]
        lo = g_s.iloc[-1]
        dumb_rows.append({
            "parent": parent.replace("lineage", "L"),
            "p_grp": hi["parent_group"],
            "hi_name": hi["sublineage_L3"].replace("lineage", "L"),
            "hi_c": int(hi["n_countries"]),
            "hi_H": float(hi["shannon_H"]),
            "lo_name": lo["sublineage_L3"].replace("lineage", "L"),
            "lo_c": int(lo["n_countries"]),
            "lo_H": max(float(lo["shannon_H"]), 0.0),
            "gap": int(hi["n_countries"]) - int(lo["n_countries"]),
        })
    df_dumb = pd.DataFrame(dumb_rows).sort_values("gap", ascending=True).reset_index(drop=True)
    y = np.arange(len(df_dumb))

    for yi, (_, r) in zip(y, df_dumb.iterrows()):
        line_c = "#C44536" if r["p_grp"] == "Global" else "#2A6F82"
        ax_h.plot([r["lo_c"], r["hi_c"]], [yi, yi], color=line_c, lw=2.4, alpha=0.55, zorder=2)
        ax_h.scatter([r["lo_c"]], [yi], s=72, color="#1C595A", edgecolor="white", lw=0.8, zorder=4)
        ax_h.scatter([r["hi_c"]], [yi], s=72, color="#C44536", edgecolor="white", lw=0.8, zorder=4)
        ax_h.text(r["lo_c"] - 1.8, yi, f"{r['lo_name']} ({r['lo_c']})", ha="right", va="center", fontsize=7.4, color="#1C595A", fontweight="bold")
        ax_h.text(r["hi_c"] + 1.8, yi, f"{r['hi_name']} ({r['hi_c']})", ha="left", va="center", fontsize=7.4, color="#9B2C2C", fontweight="bold")

    ax_h.set_yticks(y)
    ax_h.set_yticklabels([f"Within {p}" for p in df_dumb["parent"]], fontsize=8.4, fontweight="bold")
    ax_h.set_xlim(-22, 114)
    ax_h.set_xlabel("Countries colonized (narrowest vs. broadest sister sublineage)")
    ax_h.set_title("h   Within-clade nested specialization span (L3 sister divergence)",
                   loc="left", fontweight="bold", fontsize=10.2, pad=8)
    ax_h.grid(axis="x", linestyle=":", alpha=0.35)
    h_leg = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#C44536", markersize=7.5, label="Broadest L3 sister"),
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor="#1C595A", markersize=7.5, label="Most localized L3 sister"),
    ]
    ax_h.legend(handles=h_leg, loc="lower right", frameon=True, facecolor="#FAF8F4", edgecolor="#D6D3CD", fontsize=7.8)

    save_v2(fig, "Fig3_ThirdLevel_Nested_Specialization_V2")


if __name__ == "__main__":
    main()
