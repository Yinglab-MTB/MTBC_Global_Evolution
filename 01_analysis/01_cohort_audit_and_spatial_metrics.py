# -*- coding: utf-8 -*-
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import openpyxl
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, SUB14, GLOBAL_SET, ROOTED_SET, ensure_output_dirs


def build_cohort_funnel_and_breakdown():
    ensure_output_dirs()

    xlsx_fig1 = DATA_DIR / "MTBC_global_cohort_metadata_and_summaries.xlsx"
    wb1 = openpyxl.load_workbook(xlsx_fig1, read_only=True, data_only=True)

    ws1 = wb1["Sheet1"]
    n_total = 0
    years = []
    lin_cnt = Counter()
    dr_cnt = Counter()
    for i, r in enumerate(ws1.iter_rows(values_only=True)):
        if i == 0:
            continue
        if r[0] is None:
            continue
        n_total += 1
        if r[1] is not None:
            try:
                y = int(r[1])
                if 1850 <= y <= 2030:
                    years.append(y)
            except (TypeError, ValueError):
                pass
        lin_cnt[str(r[2]).strip()] += 1
        dr_cnt[str(r[3]).strip()] += 1

    wb1.close()

    df_csv = pd.read_csv(DATA_DIR / "MTBC_country_lineage_counts_coordinates.csv")
    n_geo_mapped = int(df_csv["Size"].sum())
    n_countries = len(df_csv)
    l4_geo = int(df_csv["lineage4"].sum())
    l2_geo = int(df_csv["lineage2"].sum())
    l1_geo = int(df_csv["lineage1"].sum())
    l3_geo = int(df_csv["lineage3"].sum())
    la1_geo = int(df_csv["La1"].sum())
    other_geo = n_geo_mapped - (l4_geo + l2_geo + l1_geo + l3_geo + la1_geo)

    df_bur = pd.read_csv(
        DATA_DIR / "sublineage_drug_resistance_burden.tsv", sep="\t"
    )
    n_14sub = int(df_bur["n_samples"].sum())
    n_global = int(df_bur[df_bur["sublineage"].isin(GLOBAL_SET)]["n_samples"].sum())
    n_rooted = int(df_bur[df_bur["sublineage"].isin(ROOTED_SET)]["n_samples"].sum())
    n_inter = n_14sub - n_global - n_rooted

    beast_ntax_raw = 0
    for sn in SUB14:
        p = DATA_DIR / "phylogeographic_trees" / f"{sn}_location_tree_with_trait.tree"
        txt = p.read_text(encoding="utf-8", errors="ignore")
        m = re.search(r"ntax\s*=\s*(\d+)", txt, re.I)
        if m:
            beast_ntax_raw += int(m.group(1))
    beast_ntax_postqc = beast_ntax_raw - 1

    yr_min = min(years) if years else 1900
    yr_max = max(years) if years else 2025

    bd_rows = []
    for k, v in lin_cnt.most_common():
        bd_rows.append({
            "Category": "Main_Lineage",
            "Label": k,
            "Count": v,
            "Percentage": round(100.0 * v / n_total, 3),
        })
    for k, v in dr_cnt.most_common():
        bd_rows.append({
            "Category": "Drug_Resistance_Type",
            "Label": k,
            "Count": v,
            "Percentage": round(100.0 * v / n_total, 3),
        })
    df_bd = pd.DataFrame(bd_rows)
    df_bd.to_csv(TABLES_DIR / "Table_S1b_MainLineage_and_DR_Breakdown.tsv", sep="\t", index=False)
    return df_bd


def compute_spatial_metrics():
    ensure_output_dirs()

    wb2 = openpyxl.load_workbook(
        DATA_DIR / "MTBC_sublineage_country_frequencies_levels_2_3.xlsx", read_only=True, data_only=True
    )

    def _summarize_sheet(ws):
        rows = list(ws.iter_rows(values_only=True))
        recs = [
            (str(r[0]), float(r[1]))
            for r in rows[1:]
            if r[0] is not None and r[1] is not None and isinstance(r[1], (int, float)) and r[1] > 0
        ]
        if not recs:
            return None
        total = sum(x[1] for x in recs)
        recs_s = sorted(recs, key=lambda x: x[1], reverse=True)
        ps = np.array([x[1] / total for x in recs], dtype=float)
        H = float(-np.sum(ps * np.log(ps)))
        top3_share = sum(x[1] for x in recs_s[:3]) / total
        return {
            "N_geo": int(round(total)),
            "n_countries": len(recs),
            "shannon_H": round(H, 4),
            "top1_country": f"{recs_s[0][0]} ({recs_s[0][1] / total * 100:.1f}%)",
            "top3_countries": "; ".join(f"{c} ({n / total * 100:.1f}%)" for c, n in recs_s[:3]),
            "top3_share": round(top3_share, 4),
        }

    l2_rows = []
    for sn in SUB14:
        m = _summarize_sheet(wb2[sn])
        if m is None:
            continue
        grp = "Global" if sn in GLOBAL_SET else ("Rooted" if sn in ROOTED_SET else "Intermediate")
        l2_rows.append({"sublineage": sn, "ecological_group": grp, **m})

    l3_rows = []
    for sn in wb2.sheetnames:
        if sn.startswith("lineage") and sn.count(".") == 2:
            m = _summarize_sheet(wb2[sn])
            if m is None:
                continue
            parent = ".".join(sn.split(".")[:2])
            pgrp = "Global" if parent in GLOBAL_SET else ("Rooted" if parent in ROOTED_SET else "Intermediate")
            l3_rows.append({
                "sublineage_L3": sn,
                "parent_L2": parent,
                "parent_group": pgrp,
                "N_geo": m["N_geo"],
                "n_countries": m["n_countries"],
                "shannon_H": m["shannon_H"],
                "top1_country": m["top1_country"],
                "top3_share": m["top3_share"],
            })
    wb2.close()

    df_l2 = pd.DataFrame(l2_rows)
    df_l3 = pd.DataFrame(l3_rows).sort_values(["parent_L2", "N_geo"], ascending=[True, False])

    df_l2.to_csv(TABLES_DIR / "Table_S2_spatial_metrics_L2.tsv", sep="\t", index=False)
    df_l3.to_csv(TABLES_DIR / "Table_S2c_ThirdLevel_Sublineage_Geography.tsv", sep="\t", index=False)
    return df_l2, df_l3


def main():
    build_cohort_funnel_and_breakdown()
    compute_spatial_metrics()


if __name__ == "__main__":
    main()
