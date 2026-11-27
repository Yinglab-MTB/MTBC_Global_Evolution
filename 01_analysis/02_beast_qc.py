# -*- coding: utf-8 -*-
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, SUB14, ensure_output_dirs

def _get_annot_val(annot_str, key):
    m = re.search(r"(?:^|,)" + re.escape(key) + r"=(\{[^}]+\}|\"[^\"]+\"|[^,]+)", annot_str)
    return m.group(1) if m else ""


def audit_all_14_beast_trees_and_bsps():
    ensure_output_dirs()

    tree_rows = []
    for sn in SUB14:
        tpath = DATA_DIR / "phylogeographic_trees" / f"{sn}_location_tree_with_trait.tree"
        txt = tpath.read_text(encoding="utf-8", errors="ignore")
        m_ntax = re.search(r"ntax\s*=\s*(\d+)", txt, re.I)
        ntax = int(m_ntax.group(1)) if m_ntax else 0

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

        m_tree = re.search(r"tree\s+\S+\s*=\s*(?:\[&R\]\s*)?(.+);", txt, re.S)
        tree_str = m_tree.group(1) if m_tree else ""
        annots = re.findall(r"\[&([^\]]+)\]", tree_str)
        root_a = annots[-1]

        rh = float(_get_annot_val(root_a, "height"))
        rhpd_raw = _get_annot_val(root_a, "height_95%_HPD").strip("{}").split(",")
        rhpd = f"{float(rhpd_raw[0]):.1f}-{float(rhpd_raw[1]):.1f}" if len(rhpd_raw) == 2 else ""
        lset = _get_annot_val(root_a, "location.set").strip("{}").replace('"', "").split(",")
        lprob = [float(x) for x in _get_annot_val(root_a, "location.set.prob").strip("{}").split(",")]
        pairs = sorted(zip(lset, lprob), key=lambda x: -x[1])

        heights = sorted([float(x) for x in re.findall(r"(?:\[&|,)height=([0-9.eE+-]+)", tree_str)], reverse=True)
        h2 = heights[1] if len(heights) > 1 else rh
        h_ratio = rh / h2 if h2 > 0 else 1.0

        tip_annots = re.findall(r"(?:^|[,(])[^(),\[]+\[&([^\]]+)\]", tree_str)
        tip_locs = [
            re.search(r'(?:^|,)location="?([^",\}]+)"?', a).group(1)
            for a in tip_annots
            if re.search(r"(?:^|,)location=", a)
        ]
        tc = Counter(tip_locs)

        outliers = []
        if h_ratio > 1.30:
            for match in re.finditer(r"([^(),\[:]*)\[&([^\]]+)\]:([0-9.eE+-]+)", tree_str):
                lbl, annot, blen = match.group(1), match.group(2), float(match.group(3))
                if lbl and blen > 0.50 * rh:
                    mloc = re.search(r'(?:^|,)location="?([^",\}]+)"?', annot)
                    outliers.append(f"{tmap.get(lbl, lbl)} ({mloc.group(1) if mloc else '?'}, len={blen:.1f}y)")

        top1_loc, top1_pp = pairs[0]
        if h_ratio >= 1.8:
            qc_flag = "CRITICAL_OUTLIER_RERUN_REQUIRED"
        elif h_ratio >= 1.35:
            qc_flag = "CHECK_DEEP_BASAL_BRANCHES"
        elif top1_pp >= 0.65:
            qc_flag = "HIGH_CONFIDENCE_ORIGIN"
        else:
            qc_flag = "UNCERTAIN_DEEP_ORIGIN_RAPID_RADIATION"

        flagged_str = "; ".join(outliers) if outliers else "None"
        tree_rows.append({
            "sublineage": sn,
            "BEAST_ntax": ntax,
            "root_height_yr": round(rh, 1),
            "root_height_95HPD": rhpd,
            "second_deepest_H_yr": round(h2, 1),
            "H1_H2_ratio": round(h_ratio, 2),
            "root_top1_loc": top1_loc,
            "root_top1_PP": round(top1_pp, 3),
            "root_top3_summary": "; ".join(f"{k}:{v:.3f}" for k, v in pairs[:3]),
            "tip_region_balance": "; ".join(f"{k}={v}" for k, v in tc.most_common(5)),
            "QC_status": qc_flag,
            "flagged_long_tips": flagged_str,
        })

    df_tree = pd.DataFrame(tree_rows)
    df_tree.to_csv(TABLES_DIR / "Table_S4_BEAST_Tree_QC_and_Outlier_Audit.tsv", sep="\t", index=False)

    bsp_rows = []
    for sn in SUB14:
        short = sn.replace("lineage", "L")
        p = DATA_DIR / "bayesian_skyline" / f"{short}.BSP.txt"
        lines = [l.strip() for l in p.read_text(encoding="utf-8", errors="ignore").splitlines() if l.strip()]
        h_idx = next(i for i, l in enumerate(lines) if l.startswith("Time"))
        header = lines[h_idx].split("\t")
        rows = [[float(x) if x not in ("", "NaN", "nan") else np.nan for x in l.split("\t")] for l in lines[h_idx + 1:]]
        d = pd.DataFrame(rows, columns=header).dropna(subset=["Median"]).sort_values("Time").reset_index(drop=True)

        t0, t1 = float(d["Time"].iloc[0]), float(d["Time"].iloc[-1])
        ne0, ne1 = float(d["Median"].iloc[0]), float(d["Median"].iloc[-1])
        pk_i = d["Median"].idxmax()
        t_pk, ne_pk = float(d.loc[pk_i, "Time"]), float(d.loc[pk_i, "Median"])
        ne_1950 = float(d[d["Time"] >= 1950]["Median"].iloc[0])
        ne_2000 = float(d[d["Time"] >= 2000]["Median"].iloc[0])

        bsp_rows.append({
            "sublineage": sn,
            "BSP_t_start": round(t0, 1),
            "BSP_Ne_start": round(ne0, 1),
            "BSP_t_peak": round(t_pk, 1),
            "BSP_Ne_peak": round(ne_pk, 1),
            "BSP_Ne_1950": round(ne_1950, 1),
            "BSP_Ne_2000": round(ne_2000, 1),
            "BSP_Ne_2025": round(ne1, 1),
            "BSP_fold_expansion": round(ne_pk / ne0, 1),
            "retain_vs_peak": round(ne1 / ne_pk, 4),
            "retain_vs_1950": round(ne1 / ne_1950, 4),
            "retain_vs_2000": round(ne1 / ne_2000, 4),
        })

    df_bsp = pd.DataFrame(bsp_rows)
    df_bsp.to_csv(TABLES_DIR / "Table_S2_bsp_metrics.tsv", sep="\t", index=False)
    return df_tree, df_bsp

def main():
    return audit_all_14_beast_trees_and_bsps()


if __name__ == "__main__":
    main()
