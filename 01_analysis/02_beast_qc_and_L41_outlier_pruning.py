# -*- coding: utf-8 -*-
"""
Step 02: BEAST MCC Tree Quality-Control Audit, Post-QC L4.1 Outlier Pruning
         (Isolate 2405024009), Crown-Aligned BSP Extraction & FASTA Filter CLI.

Scope:
  - Performs systematic H1/H2 tree-balance QC across all 14 BEAST MCC trees.
  - Excludes the single anomalous basal outlier isolate `2405024009` ONLY from
    `L4.1` BEAST phylodynamic (BSP) and phylogeographic (MCC tree) analyses,
    producing:
      * results/cleaned_L41_beast/lineage4.1_location_tree_with_trait_clean499.tree
      * results/cleaned_L41_beast/L4.1_clean499.BSP.txt
      * results/tables/Table_S4_BEAST_Tree_QC_and_Outlier_Audit.tsv
      * results/tables/Table_S2_bsp_metrics.tsv (for Step 03 integration)
  - Also exposes `--filter-fasta` and `--filter-meta` CLI options to exclude
    `2405024009` from raw alignment/trait tables prior to a fresh BEAST MCMC rerun.
"""

import argparse
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, CLEAN_BEAST_DIR, SUB14, ensure_output_dirs

FLAGGED_OUTLIER_IDS = {"2405024009"}


def _get_annot_val(annot_str, key):
    m = re.search(r"(?:^|,)" + re.escape(key) + r"=(\{[^}]+\}|\"[^\"]+\"|[^,]+)", annot_str)
    return m.group(1) if m else ""


def prune_L41_outlier_tree():
    """
    Prune isolated basal outlier '2405024009' from the L4.1 BEAST MCC NEXUS tree
    and promote its 499-tip sister crown node to the post-QC root.
    """
    src = DATA_DIR / "Fig. 5" / "lineage4.1_location_tree_with_trait.tree"
    dst = CLEAN_BEAST_DIR / "lineage4.1_location_tree_with_trait_clean499.tree"
    txt = src.read_text(encoding="utf-8", errors="ignore")

    txt_new = re.sub(r"ntax\s*=\s*500", "ntax=499", txt, flags=re.I)
    txt_new = re.sub(r"^\s*'?2405024009'?\s*\r?\n", "", txt_new, flags=re.M)
    txt_new = re.sub(r"^\s*1\s+'?2405024009'?,?\s*\r?\n", "", txt_new, flags=re.M)

    m_tree = re.search(r"(tree\s+\S+\s*=\s*(?:\[&R\]\s*)?)(.+)(;\s*End;)", txt_new, re.S | re.I)
    if not m_tree:
        m_tree = re.search(r"(tree\s+\S+\s*=\s*(?:\[&R\]\s*)?)(.+)(;)", txt_new, re.S | re.I)
    prefix, newick, suffix = m_tree.group(1), m_tree.group(2).strip(), m_tree.group(3)

    # Split outermost root children: (outlier_branch, crown_499_subtree)[&root_annot]
    depth = 0
    in_br = False
    split_pos = None
    outer_end = None
    for i, ch in enumerate(newick):
        if ch == "[":
            in_br = True
        elif ch == "]":
            in_br = False
        elif not in_br:
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0:
                    outer_end = i
                    break
            elif ch == "," and depth == 1:
                split_pos = i

    child1 = newick[1:split_pos].strip()
    child2 = newick[split_pos + 1:outer_end].strip()
    crown_sub = child2 if (child1.startswith("1[") or child1.startswith("1:")) else child1
    crown_clean = re.sub(r":[0-9.eE+-]+$", "", crown_sub.strip())

    txt_final = txt_new[:m_tree.start()] + prefix + crown_clean + suffix + txt_new[m_tree.end():]
    dst.write_text(txt_final, encoding="utf-8")

    # Parse crown root height for BSP alignment
    annots = re.findall(r"\[&([^\]]+)\]", crown_clean)
    rh = float(_get_annot_val(annots[-1], "height"))
    print(f"  [L4.1 Tree Pruned] Saved {dst.name} (ntax=499, crown tMRCA={rh:.1f} y)")
    return dst, rh


def align_L41_bsp_to_crown(crown_height_yr):
    """Align L4.1 Tracer BSP table to the 499-tip crown onset (2025 - crown_height)."""
    src = DATA_DIR / "Fig. 4" / "L4.1.BSP.txt"
    dst = CLEAN_BEAST_DIR / "L4.1_clean499.BSP.txt"

    raw_lines = src.read_text(encoding="utf-8", errors="ignore").splitlines()
    h_idx = next(i for i, l in enumerate(raw_lines) if l.strip().startswith("Time"))
    header = raw_lines[h_idx].strip().split("\t")
    rows = [
        [float(x) if x not in ("", "NaN", "nan") else np.nan for x in l.strip().split("\t")]
        for l in raw_lines[h_idx + 1:]
        if l.strip()
    ]
    df = pd.DataFrame(rows, columns=header).dropna(subset=["Median"]).sort_values("Time").reset_index(drop=True)

    crown_start_year = 2025.0 - crown_height_yr
    df_clean = df[df["Time"] >= crown_start_year].copy().reset_index(drop=True)

    with dst.open("w", encoding="utf-8") as fout:
        fout.write("# Post-QC L4.1 BSP (outlier 2405024009 removed; crown-aligned n=499)\n")
        fout.write("\t".join(header) + "\n")
        for _, r in df_clean.iterrows():
            fout.write("\t".join(f"{r[c]:.6g}" if pd.notna(r[c]) else "NaN" for c in header) + "\n")
    print(f"  [L4.1 BSP Aligned] Saved {dst.name} (span {df_clean['Time'].iloc[0]:.1f}-{df_clean['Time'].iloc[-1]:.1f} CE)")
    return dst


def audit_all_14_beast_trees_and_bsps():
    """Run tree-balance QC and BSP summary extraction across all 14 sublineages (using post-QC L4.1)."""
    print("[01_analysis | Step 2] Running BEAST tree QC, L4.1 pruning & BSP extraction...")
    ensure_output_dirs()

    clean_l41_tree, l41_rh = prune_L41_outlier_tree()
    clean_l41_bsp = align_L41_bsp_to_crown(l41_rh)

    # 1. Audit all 14 BEAST MCC trees
    tree_rows = []
    for sn in SUB14:
        tpath = clean_l41_tree if sn == "lineage4.1" else (DATA_DIR / "Fig. 5" / f"{sn}_location_tree_with_trait.tree")
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

        heights = sorted([float(x) for x in re.findall(r"(?:^|,)height=([0-9.eE+-]+)", tree_str)], reverse=True)
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

        flagged_str = "None (outlier 2405024009 removed)" if sn == "lineage4.1" else (
            "; ".join(outliers) if outliers else "None"
        )
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

    # 2. Extract BSP demographic trajectories & retention ratios across all 14 sublineages
    bsp_rows = []
    for sn in SUB14:
        short = sn.replace("lineage", "L")
        p = clean_l41_bsp if sn == "lineage4.1" else (DATA_DIR / "Fig. 4" / f"{short}.BSP.txt")
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
    print(f"  Saved Table_S4 ({len(df_tree)} post-QC trees; max H1/H2={df_tree['H1_H2_ratio'].max():.2f}) & BSP metrics.")
    return df_tree, df_bsp


def filter_fasta(in_path, out_path, remove_ids):
    """Remove flagged outlier taxa from a raw FASTA alignment prior to BEAST MCMC."""
    kept, removed = 0, []
    with open(in_path, "r", encoding="utf-8") as fin, open(out_path, "w", encoding="utf-8") as fout:
        write_flag = True
        for line in fin:
            if line.startswith(">"):
                header = line[1:].strip().split()[0]
                if any(rid in header for rid in remove_ids):
                    write_flag = False
                    removed.append(header)
                else:
                    write_flag = True
                    kept += 1
            if write_flag:
                fout.write(line)
    print(f"[FASTA] Kept {kept} sequences; Removed {len(removed)}: {removed}")


def filter_trait_table(in_path, out_path, remove_ids):
    """Remove flagged outlier taxa from a BEAST trait/metadata table."""
    kept, removed = 0, []
    with open(in_path, "r", encoding="utf-8") as fin, open(out_path, "w", encoding="utf-8") as fout:
        fout.write(fin.readline())
        for line in fin:
            first_col = line.strip().split("\t")[0].split(",")[0]
            if any(rid in first_col for rid in remove_ids):
                removed.append(first_col)
                continue
            fout.write(line)
            kept += 1
    print(f"[META]  Kept {kept} rows; Removed {len(removed)}: {removed}")


def main():
    parser = argparse.ArgumentParser(description="BEAST tree QC, L4.1 outlier pruning, and optional FASTA/trait filtering")
    parser.add_argument("--filter-fasta", help="Optional input FASTA alignment to filter before BEAST rerun")
    parser.add_argument("--out-fasta", help="Output cleaned FASTA alignment")
    parser.add_argument("--filter-meta", help="Optional input trait/metadata TSV or CSV to filter")
    parser.add_argument("--out-meta", help="Output cleaned trait/metadata table")
    parser.add_argument("--extra-ids", nargs="*", default=[], help="Additional sample IDs to exclude")
    args = parser.parse_args()

    if args.filter_fasta or args.filter_meta:
        remove_ids = FLAGGED_OUTLIER_IDS | set(args.extra_ids)
        if args.filter_fasta and args.out_fasta:
            filter_fasta(args.filter_fasta, args.out_fasta, remove_ids)
        if args.filter_meta and args.out_meta:
            filter_trait_table(args.filter_meta, args.out_meta, remove_ids)
    else:
        audit_all_14_beast_trees_and_bsps()


if __name__ == "__main__":
    main()
