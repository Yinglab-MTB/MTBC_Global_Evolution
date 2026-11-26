# -*- coding: utf-8 -*-
"""
Nature Editorial V2 Visual Engine:
  - Harmonized Nature Microbiology / Nature Genetics editorial color palettes
    (muted terracotta, petrol teal, ochre gold, plum, sage, warm parchment)
  - Exact mathematical Robinson world projection + curved oval boundary + graticules
  - Land-masked Ordinary Kriging with smooth contour isolines in Robinson projection
  - Raincloud (half-violin KDE + inner box + beeswarm) plotting primitive
  - Circular/Fan time-calibrated BEAST MCC tree renderer with concentric century rings
  - Pure-Python BEAST NEXUS parser (reused & enhanced)
"""

import os
import re
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from matplotlib.collections import PolyCollection, LineCollection
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Polygon, FancyBboxPatch

import sys
from pathlib import Path as _FilePath

sys.path.insert(0, str(_FilePath(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, CLEAN_BEAST_DIR, FIGURES_DIR, ASSETS_DIR, ensure_output_dirs

BASE_DIR = str(DATA_DIR)
AUDIT_DIR = str(TABLES_DIR)
V1_ASSET_DIR = str(ASSETS_DIR)
CLEAN_DIR = str(CLEAN_BEAST_DIR)
FIG_DIR = str(FIGURES_DIR)

ensure_output_dirs()

SUB14 = [
    "lineage1.1", "lineage1.2", "lineage2.1", "lineage2.2", "lineage3.1",
    "lineage4.1", "lineage4.2", "lineage4.3", "lineage4.4", "lineage4.5",
    "lineage4.6", "lineage4.7", "lineage4.8", "lineage4.9",
]
GLOBAL_SET = {"lineage2.2", "lineage4.1", "lineage4.3", "lineage4.4", "lineage4.8"}
ROOTED_SET = {"lineage2.1", "lineage4.5"}

# ==============================================================================
# NATURE EDITORIAL COLOR SYSTEM (Muted, Jewel-Toned, Colorblind-Safe)
# ==============================================================================
PAL_GROUP = {
    "Global": "#C44536",        # Deep Nature Terracotta
    "Intermediate": "#2A6F82",  # Deep Petrol Teal
    "Rooted": "#6B7280",        # Warm Slate
}
PAL_GROUP_LIGHT = {
    "Global": "#F5D0C5",
    "Intermediate": "#C5DFE6",
    "Rooted": "#DCE0E5",
}

PAL_LINEAGE = {
    "lineage4": "#C44536",   # Terracotta
    "lineage2": "#2A6F82",   # Petrol Teal
    "lineage1": "#D99B26",   # Ochre Gold
    "lineage3": "#6D597A",   # Muted Plum
    "La1":      "#3E8E7E",   # Botanical Sage
    "lineage6": "#355070",   # Deep Indigo
    "lineage5": "#84A98C",   # Soft Eucalyptus
    "La2":      "#9C6644",   # Raw Sienna
    "La3":      "#B5838D",   # Dusty Rose
    "lineage7": "#B5173B",   # Crimson Berry
    "lineage9": "#5E548E",   # Royal Violet
    "M.canetti":"#2B2D42",   # Charcoal Ink
    "lineage8": "#C8B6FF",   # Pale Lavender
    "Others":   "#A89F91",   # Warm Stone
}

PAL_REGION = {
    "East_Asia":        "#C44536",  # Terracotta
    "Southeast_Asia":   "#E59866",  # Warm Apricot
    "South_Asia":       "#D99B26",  # Ochre Gold
    "Western_Asia":     "#8E5572",  # Mauve Berry
    "Central_Asia":     "#9C6644",  # Sienna
    "Europe":           "#2E5A88",  # Classic Nature Cobalt
    "Africa":           "#2D8675",  # Deep Emerald Teal
    "Northern_America": "#5EA3B8",  # Sky Petrol
    "South_America":    "#4A4E69",  # Dusk Slate
    "Oceania":          "#8D99AE",  # Cool Pewter
    "Unknown":          "#D6D3CD",
}

PAL_AXIS = {
    "Envelope": "#C44536",
    "Cation":   "#2A6F82",
    "Kinase":   "#D99B26",
    "Adaptive": "#6D597A",
    "Tier 2":   "#6B7280",
}

# Custom Editorial Colormaps
CMAP_TERRACOTTA = LinearSegmentedColormap.from_list(
    "nature_terracotta",
    ["#F4F1EA", "#F1DEC2", "#E6AF73", "#D47244", "#A8332E", "#590D22"],
)
CMAP_TERRACOTTA.set_bad(alpha=0.0)

CMAP_BOTANICAL = LinearSegmentedColormap.from_list(
    "nature_botanical",
    ["#F3F5F2", "#D0E5DB", "#89C2B1", "#3D8C82", "#1C595A", "#0B2F36"],
)
CMAP_BOTANICAL.set_bad(alpha=0.0)

CMAP_BORDEAUX = LinearSegmentedColormap.from_list(
    "nature_bordeaux",
    ["#FAF8F4", "#F4DFB6", "#E59866", "#C44536", "#7A1C35", "#380920"],
)


def set_nature_v2_style():
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
        "axes.edgecolor": "#333840",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "xtick.color": "#333840",
        "ytick.color": "#333840",
        "text.color": "#1E2229",
        "axes.labelcolor": "#1E2229",
        "axes.labelsize": 9.2,
        "axes.titlesize": 10.2,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "legend.fontsize": 8.0,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "savefig.bbox": "tight",
    })


def save_v2(fig, stem, dpi=400):
    pdf_p = os.path.join(FIG_DIR, f"{stem}.pdf")
    png_p = os.path.join(FIG_DIR, f"{stem}.png")
    fig.savefig(pdf_p, format="pdf", bbox_inches="tight")
    fig.savefig(png_p, format="png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"  [V2] Saved {stem}.pdf & .png")


# ==============================================================================
# EXACT MATHEMATICAL ROBINSON PROJECTION (Pure Python, No Cartopy Needed!)
# ==============================================================================
_ROB_LAT = np.arange(0, 95, 5, dtype=float)
_ROB_PLEN = np.array([
    1.0000, 0.9986, 0.9954, 0.9900, 0.9822, 0.9730, 0.9600, 0.9427, 0.9216, 0.8962,
    0.8679, 0.8350, 0.7986, 0.7597, 0.7186, 0.6732, 0.6213, 0.5722, 0.5322,
])
_ROB_PDFE = np.array([
    0.0000, 0.0620, 0.1240, 0.1860, 0.2480, 0.3100, 0.3720, 0.4340, 0.4958, 0.5571,
    0.6176, 0.6769, 0.7346, 0.7903, 0.8435, 0.8936, 0.9394, 0.9761, 1.0000,
])
_F_PLEN = interp1d(_ROB_LAT, _ROB_PLEN, kind="cubic", fill_value="extrapolate")
_F_PDFE = interp1d(_ROB_LAT, _ROB_PDFE, kind="cubic", fill_value="extrapolate")


def robinson_project(lon, lat):
    """Project (lon, lat) in degrees to Robinson (x, y) coordinates."""
    lon_a = np.asarray(lon, dtype=float)
    lat_a = np.asarray(lat, dtype=float)
    abs_lat = np.clip(np.abs(lat_a), 0.0, 90.0)
    x = 0.8487 * np.deg2rad(lon_a) * _F_PLEN(abs_lat)
    y = 1.3523 * np.sign(lat_a) * _F_PDFE(abs_lat)
    return x, y


def robinson_boundary_polygon(lat_min=-58.0, lat_max=84.0):
    lats_r = np.linspace(lat_min, lat_max, 90)
    lats_l = np.linspace(lat_max, lat_min, 90)
    x_r, y_r = robinson_project(np.full_like(lats_r, 180.0), lats_r)
    x_l, y_l = robinson_project(np.full_like(lats_l, -180.0), lats_l)
    lons_t = np.linspace(180.0, -180.0, 60)
    x_t, y_t = robinson_project(lons_t, np.full_like(lons_t, lat_max))
    lons_b = np.linspace(-180.0, 180.0, 60)
    x_b, y_b = robinson_project(lons_b, np.full_like(lons_b, lat_min))
    xs = np.concatenate([x_r, x_t, x_l, x_b])
    ys = np.concatenate([y_r, y_t, y_l, y_b])
    return np.column_stack([xs, ys])


_WORLD_LONLAT_CACHE = None
_WORLD_ROB_CACHE = None


def load_world_polygons_lonlat():
    global _WORLD_LONLAT_CACHE
    if _WORLD_LONLAT_CACHE is not None:
        return _WORLD_LONLAT_CACHE
    df = pd.read_csv(os.path.join(V1_ASSET_DIR, "world_map_polygons.csv"))
    df = df[df["lat"] > -58.0].copy()
    polys = []
    for _, g in df.groupby("group"):
        if len(g) >= 3:
            pts = g[["long", "lat"]].values
            if np.max(np.abs(pts[:, 0])) <= 180.5:
                polys.append(pts)
    _WORLD_LONLAT_CACHE = polys
    return polys


def load_world_polygons_robinson():
    global _WORLD_ROB_CACHE
    if _WORLD_ROB_CACHE is not None:
        return _WORLD_ROB_CACHE
    polys = []
    for pts in load_world_polygons_lonlat():
        rx, ry = robinson_project(pts[:, 0], pts[:, 1])
        polys.append(np.column_stack([rx, ry]))
    _WORLD_ROB_CACHE = polys
    return polys


def draw_robinson_basemap(ax, land_fc="#ECEAE3", land_ec="#9C988E", land_lw=0.35,
                          ocean_fc="#F4F8FA", graticule=True, border_ec="#B8C2CC"):
    """Draw a sleek Nature-editorial Robinson-projected world map."""
    bnd = robinson_boundary_polygon(-58.0, 84.0)
    ocean = Polygon(bnd, closed=True, facecolor=ocean_fc, edgecolor=border_ec, linewidth=0.7, zorder=0)
    ax.add_patch(ocean)

    if graticule:
        for lat_g in [-30, 0, 30, 60]:
            ls = np.linspace(-180, 180, 120)
            gx, gy = robinson_project(ls, np.full_like(ls, lat_g))
            ax.plot(gx, gy, color="#DCE5EC", linewidth=0.45, linestyle="-" if lat_g == 0 else ":", zorder=0.5)
        for lon_g in [-120, -60, 0, 60, 120]:
            la = np.linspace(-58, 84, 80)
            gx, gy = robinson_project(np.full_like(la, lon_g), la)
            ax.plot(gx, gy, color="#DCE5EC", linewidth=0.4, linestyle=":", zorder=0.5)

    polys = load_world_polygons_robinson()
    pc = PolyCollection(polys, facecolor=land_fc, edgecolor=land_ec, linewidth=land_lw, zorder=1)
    ax.add_collection(pc)

    ax.set_xlim(bnd[:, 0].min() - 0.06, bnd[:, 0].max() + 0.06)
    ax.set_ylim(bnd[:, 1].min() - 0.04, bnd[:, 1].max() + 0.04)
    ax.set_aspect("equal")
    ax.axis("off")
    return ocean


_ROB_MASK_CACHE = {}


def get_robinson_land_grid(nx=220, ny=105):
    """Precompute lon/lat grid mapped onto Robinson projection + exact land mask."""
    key = (nx, ny)
    if key in _ROB_MASK_CACHE:
        return _ROB_MASK_CACHE[key]
    lons = np.linspace(-178, 179, nx)
    lats = np.linspace(-56, 82, ny)
    LL, AA = np.meshgrid(lons, lats)
    pts = np.column_stack([LL.ravel(), AA.ravel()])
    mask = np.zeros(len(pts), dtype=bool)
    for poly in load_world_polygons_lonlat():
        xmin, ymin = poly.min(axis=0)
        xmax, ymax = poly.max(axis=0)
        if (xmax - xmin) * (ymax - ymin) < 0.7:
            continue
        in_bb = (pts[:, 0] >= xmin) & (pts[:, 0] <= xmax) & (pts[:, 1] >= ymin) & (pts[:, 1] <= ymax)
        if np.any(in_bb):
            mask[in_bb] |= Path(poly).contains_points(pts[in_bb])
    mask_2d = mask.reshape(ny, nx)
    RX, RY = robinson_project(LL, AA)
    _ROB_MASK_CACHE[key] = (LL, AA, RX, RY, mask_2d)
    return LL, AA, RX, RY, mask_2d


def kriging_robinson_surface(obs_lon, obs_lat, obs_val, nx=210, ny=100, range_deg=24.0, nugget=0.06):
    """Ordinary Kriging evaluated on land pixels and returned on Robinson (RX, RY) mesh."""
    LL, AA, RX, RY, mask = get_robinson_land_grid(nx, ny)
    coords = np.column_stack([obs_lon, obs_lat])
    vals = np.asarray(obs_val, dtype=float)
    n = len(vals)
    D = cdist(coords, coords)
    C = (1.0 - nugget) * np.exp(-D / range_deg) + nugget * np.eye(n)
    K = np.ones((n + 1, n + 1))
    K[:n, :n] = C
    K[n, n] = 0.0
    try:
        K_inv = np.linalg.inv(K)
    except np.linalg.LinAlgError:
        K_inv = np.linalg.pinv(K)

    land_pts = np.column_stack([LL[mask], AA[mask]])
    D_pred = cdist(land_pts, coords)
    c_pred = (1.0 - nugget) * np.exp(-D_pred / range_deg)
    rhs = np.hstack([c_pred, np.ones((len(land_pts), 1))])
    pred = np.clip((rhs @ K_inv[:, :n]) @ vals, 0.0, max(vals.max(), 1e-6))

    grid = np.full((ny, nx), np.nan)
    grid[mask] = pred
    return RX, RY, grid


# ==============================================================================
# RAINCLOUD PLOT PRIMITIVE (Half-Violin KDE + Sleek Inner Box + Jitter Dots)
# ==============================================================================
def draw_raincloud(ax, data_groups, positions, colors, dot_color_lists=None, width=0.34):
    rng = np.random.default_rng(2026)
    for idx, (vals, pos, col) in enumerate(zip(data_groups, positions, colors)):
        vals = np.asarray(vals, dtype=float)
        # 1. Right-side half-violin KDE
        if len(vals) >= 3 and np.std(vals) > 1e-6:
            kde = gaussian_kde(vals, bw_method=0.55)
            lo = max(vals.min() - 0.18 * (vals.max() - vals.min() + 1), 0)
            hi = vals.max() + 0.18 * (vals.max() - vals.min() + 1)
            ys = np.linspace(lo, hi, 120)
            dens = kde(ys)
            dens = dens / dens.max() * width
            ax.fill_betweenx(ys, pos + 0.04, pos + 0.04 + dens, facecolor=col, alpha=0.28, edgecolor="none", zorder=2)
            ax.plot(pos + 0.04 + dens, ys, color=col, linewidth=1.15, alpha=0.85, zorder=3)

        # 2. Compact sleek boxplot centered slightly left
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        iqr = q3 - q1
        w_lo = vals[vals >= q1 - 1.5 * iqr].min()
        w_hi = vals[vals <= q3 + 1.5 * iqr].max()
        bx = pos - 0.05
        ax.plot([bx, bx], [w_lo, w_hi], color="#2D3748", linewidth=1.1, zorder=3)
        rect = FancyBboxPatch((bx - 0.055, q1), 0.11, max(q3 - q1, 0.4),
                              boxstyle="round,pad=0.0,rounding_size=0.02",
                              facecolor=col, alpha=0.55, edgecolor="#1E2229", linewidth=0.9, zorder=4)
        ax.add_patch(rect)
        ax.plot([bx - 0.055, bx + 0.055], [med, med], color="#111111", linewidth=2.0, zorder=5)

        # 3. Individual data points on left rain strip
        dcols = dot_color_lists[idx] if dot_color_lists is not None else [col] * len(vals)
        xs = pos - 0.19 + rng.normal(0, 0.032, len(vals))
        for x_pt, y_pt, dc in zip(xs, vals, dcols):
            ax.scatter(x_pt, y_pt, s=46, color=dc, edgecolor="white", linewidth=0.75, zorder=6)


# ==============================================================================
# BEAST MCC TREE PARSER & CIRCULAR FAN TREE RENDERER
# ==============================================================================
def parse_beast_nexus_tree(tree_path):
    txt = open(tree_path, encoding="utf-8", errors="ignore").read()
    tmap = {}
    mt = re.search(r"Translate\s+(.*?);", txt, re.S | re.I)
    if mt:
        for line in mt.group(1).splitlines():
            line = line.strip().rstrip(",")
            if not line:
                continue
            parts = line.split(None, 1)
            if len(parts) == 2:
                tmap[parts[0]] = parts[1].strip("'\"")

    m = re.search(r"tree\s+\S+\s*=\s*(?:\[&R\]\s*)?(.+);", txt, re.S)
    s = m.group(1).strip()

    nodes = {}
    nid_c = [0]

    def new_node(parent=None):
        nid = nid_c[0]
        nid_c[0] += 1
        nodes[nid] = {
            "id": nid, "parent": parent, "children": [], "length": 0.0,
            "height": 0.0, "location": "Unknown", "loc_prob": 1.0,
            "is_tip": False, "label": "",
        }
        return nid

    def parse_annot(annot_str):
        loc, prob, h = "Unknown", 1.0, None
        ml = re.search(r'(?:^|,)location="?([^",\}]+)"?', annot_str)
        if ml:
            loc = ml.group(1)
        mp = re.search(r"(?:^|,)location\.prob=([0-9.eE+-]+)", annot_str)
        if mp:
            prob = float(mp.group(1))
        mh = re.search(r"(?:^|,)height=([0-9.eE+-]+)", annot_str)
        if mh:
            h = float(mh.group(1))
        return loc, prob, h

    i, n = 0, len(s)
    root = new_node(None)
    curr = root
    stack = []
    while i < n:
        ch = s[i]
        if ch == "(":
            child = new_node(curr)
            nodes[curr]["children"].append(child)
            stack.append(curr)
            curr = child
            i += 1
        elif ch == ",":
            parent = stack[-1]
            child = new_node(parent)
            nodes[parent]["children"].append(child)
            curr = child
            i += 1
        elif ch == ")":
            curr = stack.pop()
            i += 1
        elif ch == "[":
            j = s.find("]", i)
            annot = s[i + 2 : j] if s[i + 1 : i + 2] == "&" else s[i + 1 : j]
            loc, prob, h = parse_annot(annot)
            nodes[curr]["location"] = loc
            nodes[curr]["loc_prob"] = prob
            if h is not None:
                nodes[curr]["height"] = h
            i = j + 1
        elif ch == ":":
            i += 1
            if i < n and s[i] == "[":
                i = s.find("]", i) + 1
            j = i
            while j < n and s[j] not in ",()[] ;":
                j += 1
            try:
                nodes[curr]["length"] = float(s[i:j])
            except Exception:
                pass
            i = j
        elif ch in " \t\r\n;":
            i += 1
        else:
            j = i
            while j < n and s[j] not in ",():[]":
                j += 1
            lbl = s[i:j].strip().strip("'\"")
            nodes[curr]["label"] = tmap.get(lbl, lbl)
            nodes[curr]["is_tip"] = True
            i = j

    if len(nodes[root]["children"]) == 1:
        real_root = nodes[root]["children"][0]
        nodes[real_root]["parent"] = None
        root = real_root

    def count_tips(u):
        if not nodes[u]["children"]:
            nodes[u]["n_tips"] = 1
            return 1
        c = sum(count_tips(v) for v in nodes[u]["children"])
        nodes[u]["children"].sort(key=lambda v: nodes[v]["n_tips"])
        nodes[u]["n_tips"] = c
        return c

    count_tips(root)
    yc = [0]

    def assign_y(u):
        if not nodes[u]["children"]:
            nodes[u]["y"] = float(yc[0])
            yc[0] += 1
        else:
            for v in nodes[u]["children"]:
                assign_y(v)
            ys = [nodes[v]["y"] for v in nodes[u]["children"]]
            nodes[u]["y"] = 0.5 * (min(ys) + max(ys))

    assign_y(root)
    max_h = max(nd["height"] for nd in nodes.values())
    return nodes, root, max_h


def draw_circular_beast_tree(ax, tree_path, title_text, badge_text, accent_col, open_gap_deg=26.0):
    """Draw a sleek Nature-style circular/fan time-calibrated BEAST MCC tree."""
    nodes, root, max_h = parse_beast_nexus_tree(tree_path)
    ntips = nodes[root]["n_tips"]

    # Map y in [0, ntips-1] -> angle in radians leaving a clean gap at top-right for time labels
    start_rad = np.deg2rad(open_gap_deg / 2.0)
    span_rad = np.deg2rad(360.0 - open_gap_deg)
    r_inner = 0.14 * max_h  # sleek central ring hole

    def node_polar(nd):
        r = r_inner + (max_h - nd["height"])
        theta = start_rad + (nd["y"] / max(ntips - 1, 1)) * span_rad
        return r, theta

    # Concentric subtle time rings (e.g., 25%, 50%, 75%, 100% of max_h)
    thetas_ring = np.linspace(start_rad, start_rad + span_rad, 160)
    for frac in [0.25, 0.5, 0.75, 1.0]:
        rr = r_inner + frac * max_h
        ax.plot(rr * np.cos(thetas_ring), rr * np.sin(thetas_ring),
                color="#E5E7EB", linewidth=0.55, linestyle=":" if frac < 1.0 else "-", zorder=0)

    segs, cols = [], []
    for nid, nd in nodes.items():
        pid = nd["parent"]
        if pid is None:
            continue
        pnd = nodes[pid]
        r_c, th_c = node_polar(nd)
        r_p, th_p = node_polar(pnd)
        # Radial branch (parent radius -> child radius at child angle)
        segs.append([(r_p * np.cos(th_c), r_p * np.sin(th_c)),
                     (r_c * np.cos(th_c), r_c * np.sin(th_c))])
        cols.append(PAL_REGION.get(nd["location"], "#999999"))
        # Arc connector at parent radius from th_p to th_c
        n_arc = max(3, int(abs(th_c - th_p) / 0.04) + 2)
        th_arc = np.linspace(th_p, th_c, n_arc)
        arc_pts = np.column_stack([r_p * np.cos(th_arc), r_p * np.sin(th_arc)])
        for k in range(len(arc_pts) - 1):
            segs.append([arc_pts[k], arc_pts[k + 1]])
            cols.append(PAL_REGION.get(pnd["location"], "#999999"))

    lc = LineCollection(segs, colors=cols, linewidths=0.65, zorder=2)
    ax.add_collection(lc)

    # Tip outer color ring (makes regional dominance jump out instantly!)
    tip_r = r_inner + max_h * 1.045
    tip_xs, tip_ys, tip_cs = [], [], []
    for nd in nodes.values():
        if nd["is_tip"]:
            _, th = node_polar(nd)
            tip_xs.append(tip_r * np.cos(th))
            tip_ys.append(tip_r * np.sin(th))
            tip_cs.append(PAL_REGION.get(nd["location"], "#999999"))
    ax.scatter(tip_xs, tip_ys, s=3.8, c=tip_cs, linewidths=0, zorder=3)

    # Central root dot
    r0, th0 = node_polar(nodes[root])
    root_col = PAL_REGION.get(nodes[root]["location"], "#333333")
    ax.scatter([r0 * np.cos(th0)], [r0 * np.sin(th0)], s=26, color=root_col, edgecolor="black", linewidth=0.7, zorder=5)

    # Subtle tMRCA scale label inside open gap
    rmax = r_inner + max_h * 1.13
    ax.text(r_inner + max_h * 0.52, 0.0, f"{int(round(max_h))} y", ha="center", va="center",
            fontsize=6.6, color="#6B7280", fontweight="bold")
    ax.set_xlim(-rmax, rmax)
    ax.set_ylim(-rmax, rmax)
    ax.set_aspect("equal")
    ax.axis("off")

    ax.set_title(title_text, fontsize=9.0, fontweight="bold", color=accent_col, pad=3)
    ax.text(0.5, -0.04, badge_text, transform=ax.transAxes, ha="center", va="top", fontsize=7.4,
            color="#1E2229", bbox=dict(boxstyle="round,pad=0.22", facecolor="#FAF8F4", edgecolor=accent_col, lw=0.8))
