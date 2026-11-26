# -*- coding: utf-8 -*-
"""
Fig. 2 (Nature Editorial V2): Robinson-projected Ordinary Kriging contour maps
for all 14 major second-level MTBC sublineages + integrated ecological breadth strip.
"""

import os
import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
from matplotlib.patches import Patch

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from nature_v2_style import (
    BASE_DIR, AUDIT_DIR, PAL_GROUP, CMAP_TERRACOTTA,
    set_nature_v2_style, save_v2, robinson_project, draw_robinson_basemap, kriging_robinson_surface,
)


def main():
    print("[Fig. 2 V2] Generating Robinson-projected Kriging contour maps (14 sublineages)...")
    set_nature_v2_style()

    df_s2 = pd.read_csv(
        os.path.join(AUDIT_DIR, "Table_S2_Sublineage_Geography_BSP_DR_Integrated.tsv"), sep="\t"
    ).set_index("sublineage")
    df_c = pd.read_csv(os.path.join(BASE_DIR, "Fig. 1", "data.csv"))
    coord_map = dict(zip(df_c["Country"], zip(df_c["Longitude"], df_c["Latitude"])))

    wb = openpyxl.load_workbook(
        os.path.join(BASE_DIR, "Fig. 2", "二级谱系分布频率.xlsx"), read_only=True, data_only=True
    )

    ordered_14 = [
        "lineage2.2", "lineage4.1", "lineage4.3", "lineage4.4", "lineage4.8",
        "lineage2.1", "lineage4.5",
        "lineage1.1", "lineage1.2", "lineage3.1",
        "lineage4.2", "lineage4.6", "lineage4.7", "lineage4.9",
    ]

    fig = plt.figure(figsize=(16.8, 12.6))
    gs = gridspec.GridSpec(5, 3, hspace=0.22, wspace=0.06)

    for idx, sn in enumerate(ordered_14):
        ax = fig.add_subplot(gs[idx // 3, idx % 3])
        ocean = draw_robinson_basemap(ax, land_fc="#F2EFE8", land_ec="#A39E93", land_lw=0.26,
                                      ocean_fc="#F4F8FA", graticule=False, border_ec="#CBD5E0")

        ws = wb[sn]
        rows = list(ws.iter_rows(values_only=True))
        obs_counts = {}
        for r in rows[1:]:
            if r[0] is not None and r[1] is not None and isinstance(r[1], (int, float)):
                obs_counts[str(r[0])] = float(r[1])
        tot = sum(obs_counts.values())

        lons_o, lats_o, vals_o = [], [], []
        for ctry, (lon, lat) in coord_map.items():
            if lat < -56:
                continue
            lons_o.append(lon)
            lats_o.append(lat)
            vals_o.append(obs_counts.get(ctry, 0.0) / max(tot, 1.0))

        vmax = max(np.max(vals_o), 0.05)
        RX, RY, grid = kriging_robinson_surface(lons_o, lats_o, vals_o, nx=190, ny=92, range_deg=24.0, nugget=0.06)
        grid_norm = np.clip(grid / vmax, 0.0, 1.0)

        qm = ax.pcolormesh(RX, RY, grid_norm, cmap=CMAP_TERRACOTTA, vmin=0.0, vmax=1.0,
                           shading="gouraud", zorder=2, alpha=0.95)
        qm.set_clip_path(ocean)

        # Delicate topographic contour isolines for editorial richness
        cs = ax.contour(RX, RY, grid_norm, levels=[0.25, 0.50, 0.75], colors="#590D22",
                        linewidths=0.35, alpha=0.45, zorder=2.6)

        # Subtle coastline redraw over Kriging surface
        from matplotlib.collections import PolyCollection
        from nature_v2_style import load_world_polygons_robinson
        ax.add_collection(PolyCollection(load_world_polygons_robinson(), facecolor="none",
                                         edgecolor="#6E6A61", linewidth=0.22, zorder=3))

        # Pinpoint dots on present countries
        nz = np.array(vals_o) > 0
        px, py = robinson_project(np.array(lons_o)[nz], np.array(lats_o)[nz])
        ax.scatter(px, py, s=4.0, color="#1E2229", alpha=0.42, linewidths=0, zorder=3.5)

        meta = df_s2.loc[sn]
        grp = meta["ecological_group"]
        col = PAL_GROUP[grp]
        short = sn.replace("lineage", "L")

        # Editorial pill badge in top-left corner of map
        ax.text(0.01, 1.02, f" {short} ", transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8.3, fontweight="bold", color="white",
                bbox=dict(boxstyle="round,pad=0.18", facecolor=col, edgecolor="none"))
        ax.text(0.145, 1.02, f"{grp}  |  {int(meta['n_countries'])} countries · H = {meta['shannon_H']:.2f} · peak {vmax*100:.0f}%",
                transform=ax.transAxes, ha="left", va="bottom", fontsize=7.6, color="#2B2D42")

    wb.close()

    # 15th slot: Shared colorbar + compact ecological breadth summary strip
    ax_leg = fig.add_subplot(gs[4, 2])
    ax_leg.axis("off")

    cax = ax_leg.inset_axes([0.08, 0.78, 0.84, 0.11])
    sm = plt.cm.ScalarMappable(cmap=CMAP_TERRACOTTA, norm=Normalize(vmin=0, vmax=1))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
    cb.set_ticklabels(["0%", "25%", "50%", "75%", "100%"])
    cb.ax.tick_params(labelsize=7.6)
    cb.set_label("Relative interpolated frequency (contour lines: 25 / 50 / 75%)", fontsize=7.6, labelpad=3)

    # Mini horizontal group-comparison range strip inside legend panel
    sax = ax_leg.inset_axes([0.24, 0.08, 0.70, 0.44])
    for yi, (grp, lbl) in enumerate([("Rooted", "Rooted (n=2)"), ("Intermediate", "Regional (n=7)"), ("Global", "Global (n=5)")]):
        vals = df_s2[df_s2["ecological_group"] == grp]["n_countries"].values
        c = PAL_GROUP[grp]
        sax.plot([vals.min(), vals.max()], [yi, yi], color=c, lw=2.2, alpha=0.5, zorder=2)
        sax.scatter(vals, np.full_like(vals, yi, dtype=float), s=28, color=c, edgecolor="white", lw=0.5, zorder=3)
        sax.scatter([np.median(vals)], [yi], s=55, marker="D", color=c, edgecolor="#111111", lw=0.7, zorder=4)
    sax.set_yticks([0, 1, 2])
    sax.set_yticklabels(["Rooted (n=2)", "Regional (n=7)", "Global (n=5)"], fontsize=7.6, fontweight="bold")
    sax.set_xlim(5, 96)
    sax.set_xlabel("Countries colonized (MWU P = 0.0033)", fontsize=7.6)
    sax.tick_params(axis="x", labelsize=7.4)
    sax.grid(axis="x", linestyle=":", alpha=0.4)

    fig.suptitle("Global Kriging spatial architecture of 14 major second-level MTBC sublineages",
                 fontsize=12.5, fontweight="bold", y=0.985)
    save_v2(fig, "Fig2_SecondLevel_Sublineage_Geography_V2")


if __name__ == "__main__":
    main()
