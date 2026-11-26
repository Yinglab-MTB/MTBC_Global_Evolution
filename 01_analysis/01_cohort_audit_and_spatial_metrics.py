# -*- coding: utf-8 -*-
"""
Step 01: Four-Stage Cohort Filtering Audit & Multi-Level Spatial Niche Metrics.

Reads raw global MTBC metadata (Fig. 1) and country-level sublineage frequency
workbooks (Fig. 2 / Fig. 3) to produce:
  - Table_S1_Cohort_Filtering_Funnel.tsv
  - Table_S1b_MainLineage_and_DR_Breakdown.tsv
  - Table_S2_spatial_metrics_L2.tsv (intermediate L2 geography table for Step 03)
  - Table_S2c_ThirdLevel_Sublineage_Geography.tsv
"""

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
    """Audit the 4-stage hierarchical cohort funnel and taxonomic/WHO-DR composition."""
    print("[01_analysis | Step 1a] Auditing global cohort funnel & lineage/DR composition...")
    ensure_output_dirs()

    xlsx_fig1 = DATA_DIR / "Fig. 1" / "data.xlsx"
    wb1 = openpyxl.load_workbook(xlsx_fig1, read_only=True, data_only=True)

    # Sheet12 stores the complete 162,869-isolate curated global table
    # Columns: (target_sample, Run, Year, Country, main_lineage, sub_lineage, drtype, ...)
    ws12 = wb1["Sheet12"]
    n_total = 0
    years = []
    lin_cnt = Counter()
    dr_cnt = Counter()
    for i, r in enumerate(ws12.iter_rows(values_only=True)):
        if i == 0:
            continue
        if r[1] is None:
            continue
        n_total += 1
        if r[2] is not None:
            try:
                y = int(r[2])
                if 1850 <= y <= 2030:
                    years.append(y)
            except (TypeError, ValueError):
                pass
        lin_cnt[str(r[4]).strip()] += 1
        dr_cnt[str(r[6]).strip()] += 1

    # Stage 3a: metadata-complete isolate pool (Sheet13)
    ws13 = wb1["Sheet13"]
    n_sheet13 = sum(
        1 for i, r in enumerate(ws13.iter_rows(values_only=True))
        if i > 0 and r[0] is not None
    )
    wb1.close()

    # Stage 2: country-geolocated global atlas (120 countries)
    df_csv = pd.read_csv(DATA_DIR / "Fig. 1" / "data.csv")
    n_geo_mapped = int(df_csv["Size"].sum())
    n_countries = len(df_csv)
    l4_geo = int(df_csv["lineage4"].sum())
    l2_geo = int(df_csv["lineage2"].sum())
    l1_geo = int(df_csv["lineage1"].sum())
    l3_geo = int(df_csv["lineage3"].sum())
    la1_geo = int(df_csv["La1"].sum())
    other_geo = n_geo_mapped - (l4_geo + l2_geo + l1_geo + l3_geo + la1_geo)

    # Stage 3b: 14-sublineage comparative & drug-resistance cohort
    df_bur = pd.read_csv(
        DATA_DIR / "Extended Data Fig. 4" / "sublineage_resistance_burden.tsv", sep="\t"
    )
    n_14sub = int(df_bur["n_samples"].sum())
    n_global = int(df_bur[df_bur["sublineage"].isin(GLOBAL_SET)]["n_samples"].sum())
    n_rooted = int(df_bur[df_bur["sublineage"].isin(ROOTED_SET)]["n_samples"].sum())
    n_inter = n_14sub - n_global - n_rooted

    # Stage 4: BEAST phylodynamic cohort (post-QC removes L4.1 basal outlier 2405024009: 6,713 -> 6,712)
    beast_ntax_raw = 0
    for sn in SUB14:
        p = DATA_DIR / "Fig. 5" / f"{sn}_location_tree_with_trait.tree"
        txt = p.read_text(encoding="utf-8", errors="ignore")
        m = re.search(r"ntax\s*=\s*(\d+)", txt, re.I)
        if m:
            beast_ntax_raw += int(m.group(1))
    beast_ntax_postqc = beast_ntax_raw - 1  # L4.1: 500 -> 499

    yr_min = min(years) if years else 1900
    yr_max = max(years) if years else 2025

    funnel_rows = [
        {
            "Stage": "Stage 1: Global Curated MTBC Genomes",
            "N_Genomes": n_total,
            "Scope_Details": f"All human (L1-L9), animal (La1-La3) lineages + M. canettii; years {yr_min}-{yr_max}",
            "Source_File": "Fig. 1/data.xlsx (Sheet12)",
            "Downstream_Use": "Global genomic resource summary (Fig. 1)",
        },
        {
            "Stage": "Stage 2: Country-Geolocated Global Atlas",
            "N_Genomes": n_geo_mapped,
            "Scope_Details": (
                f"Resolved country coordinates across {n_countries} countries "
                f"(L4={l4_geo:,}; L2={l2_geo:,}; L1={l1_geo:,}; L3={l3_geo:,}; "
                f"La1={la1_geo:,}; others={other_geo:,})"
            ),
            "Source_File": "Fig. 1/data.csv & Fig. 2/3 Excel",
            "Downstream_Use": "Macro-biogeography (Fig. 1) & 2nd/3rd-level Kriging spatial mapping (Fig. 2-3)",
        },
        {
            "Stage": "Stage 3a: Metadata-Complete Sublineage Subset",
            "N_Genomes": n_sheet13,
            "Scope_Details": (
                "Non-redundant metadata-complete isolate pool with resolved Year, Country, "
                "Region, Sub-region, Sub-lineage & Drug status"
            ),
            "Source_File": "Fig. 1/data.xlsx (Sheet13)",
            "Downstream_Use": "Sampling pool for phylodynamics & comparative genomics",
        },
        {
            "Stage": "Stage 3b: 14-Sublineage Comparative & DR Cohort",
            "N_Genomes": n_14sub,
            "Scope_Details": (
                f"14 major 2nd-level sublineages: Global (5 sublineages, n={n_global}), "
                f"Rooted (L2.1+L4.5, n={n_rooted}), Intermediate (7 sublineages, n={n_inter})"
            ),
            "Source_File": "Extended Data Fig. 4 & Fig. 6",
            "Downstream_Use": "Convergent positive selection scan (Fig. 6) & DR mutation enrichment (Ext. Data Fig. 4)",
        },
        {
            "Stage": "Stage 4: Region-Balanced BEAST Phylodynamic Cohort",
            "N_Genomes": beast_ntax_postqc,
            "Scope_Details": (
                "Up to 500 geographically stratified isolates per sublineage "
                "(11 sublineages x 500 + L4.1 x 499 [after removing basal outlier 2405024009] "
                "+ L4.7 x 461 + L2.1 x 252)"
            ),
            "Source_File": "Fig. 4 (*.BSP.txt) & Fig. 5 (*.tree)",
            "Downstream_Use": "Bayesian Skyline Plots (Fig. 4) & discrete phylogeographic reconstruction (Fig. 5)",
        },
    ]
    df_funnel = pd.DataFrame(funnel_rows)
    df_funnel.to_csv(TABLES_DIR / "Table_S1_Cohort_Filtering_Funnel.tsv", sep="\t", index=False)

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
    print(f"  Saved Table_S1 (Stage 1={n_total:,}, Stage 2={n_geo_mapped:,}, Stage 4={beast_ntax_postqc:,}) & Table_S1b.")
    return df_funnel, df_bd


def compute_spatial_metrics():
    """Calculate country richness, Shannon entropy (H), and top-3 share for L2 (14) and L3 (30) sublineages."""
    print("[01_analysis | Step 1b] Computing 2nd-level (14) and 3rd-level (30) spatial niche metrics...")
    ensure_output_dirs()

    wb2 = openpyxl.load_workbook(
        DATA_DIR / "Fig. 2" / "二级谱系分布频率.xlsx", read_only=True, data_only=True
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

    # 14 second-level sublineages
    l2_rows = []
    for sn in SUB14:
        m = _summarize_sheet(wb2[sn])
        if m is None:
            continue
        grp = "Global" if sn in GLOBAL_SET else ("Rooted" if sn in ROOTED_SET else "Intermediate")
        l2_rows.append({"sublineage": sn, "ecological_group": grp, **m})

    # 30 third-level sublineages
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
    print(f"  Saved L2 spatial metrics ({len(df_l2)} clades) & Table_S2c ({len(df_l3)} L3 clades).")
    return df_l2, df_l3


def main():
    build_cohort_funnel_and_breakdown()
    compute_spatial_metrics()


if __name__ == "__main__":
    main()
