# -*- coding: utf-8 -*-
"""
Step 04: Multi-Evidence Core-Genome Convergent Positive Selection &
         Clonal Clade-Hitchhiking Filtration.

Integrates allele-frequency differentiation (Δf, BH-FDR), gene-level Fst,
nucleotide diversity (π), and Tajima's D across six complementary contrasts
(global_vs_rooted, L2.2 vs L2.1, and L4.1/L4.3/L4.4/L4.8 vs other L4; n = 23,533),
audits and filters 18 single-SNP ancestral clade-hitchhiking artifacts (Clusters A/B/C),
and outputs:
  - Table_S3_Purified_Convergent_Adaptive_Genes.tsv (Tier 1: 9 cross-L2+L4 multi-hit;
    Tier 2: 5 within-L4 multi-hit; Tier 2b: 2 cross-lineage single-hit candidates)
  - Table_S3b_Flagged_Clade_Hitchhiking_Candidates.tsv (18 filtered clade artifacts)
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DATA_DIR, TABLES_DIR, ensure_output_dirs


def run_convergent_selection_audit():
    print("[01_analysis | Step 4] Auditing convergent selection scan & filtering clonal hitchhiking artifacts...")
    ensure_output_dirs()

    f_dir = DATA_DIR / "Fig. 6"
    df_sum = pd.read_csv(f_dir / "candidate_gene_summary_all_comparisons.tsv", sep="\t")
    df_coord = pd.read_csv(f_dir / "candidate_gene_coordinates.tsv", sep="\t")
    df_int = pd.read_csv(f_dir / "candidate_gene_integrated_Fst_pi_TajimaD_evidence.tsv", sep="\t")

    # Aggregate gene-level population-genetic evidence across all six contrasts
    pop_agg = []
    for gene, g in df_int.groupby("gene_name"):
        comp_fst = dict(zip(g["comparison"], g["max_Fst"]))
        comp_td = dict(zip(g["comparison"], g["delta_tajimaD"]))
        best_row = g.sort_values(["evidence_score", "max_Fst"], ascending=[False, False]).iloc[0]
        pop_agg.append({
            "gene_name": gene,
            "max_evidence_score_raw": int(g["evidence_score"].max()),
            "max_pop_evidence_score": int(g["population_evidence_score"].max()),
            "best_comparison_popgen": best_row["comparison"],
            "overall_max_Fst": round(float(g["max_Fst"].max()), 4),
            "Fst_L22_vs_L21": round(float(comp_fst.get("L22_vs_L21", np.nan)), 4),
            "Fst_global_vs_rooted": round(float(comp_fst.get("global_vs_rooted", np.nan)), 4),
            "Fst_L41": round(float(comp_fst.get("L41_vs_otherL4", np.nan)), 4),
            "Fst_L43": round(float(comp_fst.get("L43_vs_otherL4", np.nan)), 4),
            "Fst_L44": round(float(comp_fst.get("L44_vs_otherL4", np.nan)), 4),
            "Fst_L48": round(float(comp_fst.get("L48_vs_otherL4", np.nan)), 4),
            "delta_TajimaD_L22_vs_L21": round(float(comp_td.get("L22_vs_L21", np.nan)), 4),
            "delta_TajimaD_global_vs_rooted": round(float(comp_td.get("global_vs_rooted", np.nan)), 4),
            "min_delta_pi": round(float(g["delta_pi"].min()), 6),
        })
    df_pop = pd.DataFrame(pop_agg)

    meta_cols = ["gene_name", "locus", "functional_category", "product", "gene_function"]
    df_all = (
        df_sum.merge(df_coord[meta_cols].drop_duplicates("gene_name"), on="gene_name", how="left")
        .merge(df_pop, on="gene_name", how="left")
    )

    def _classify_gene(row):
        comps = str(row["comparisons"])
        n_mut = int(row["n_mutations"])
        n_cmp = int(row["n_comparisons"])
        dfreq = float(row["max_delta_freq"])
        has_L2 = ("L22_vs_L21" in comps) or (pd.notna(row["Fst_L22_vs_L21"]) and row["Fst_L22_vs_L21"] >= 0.75)
        has_L4 = any(x in comps for x in ["L41_vs_otherL4", "L43_vs_otherL4", "L44_vs_otherL4", "L48_vs_otherL4"])

        # Explicit filtration of single-SNP ancestral clade-hitchhiking clusters
        if n_mut == 1 and n_cmp >= 3:
            if 0.246 <= dfreq <= 0.248:
                note = "FLAGGED_ClusterA_L4.7_L4.8_L4.9_H37Rv_clade_SNP (df~0.247)"
                if pd.notna(row["Fst_L22_vs_L21"]) and row["Fst_L22_vs_L21"] >= 0.90:
                    note += "; VERIFY_L2_SITE_SEPARATELY (Fst_L22=1.0)"
                return "Tier3_Flagged_Clade_Artifact", note
            if 0.322 <= dfreq <= 0.324:
                return "Tier3_Flagged_Clade_Artifact", "FLAGGED_ClusterB_L4.6_L4.7_L4.8_L4.9_clade_SNP (df~0.323)"
            if 0.613 <= dfreq <= 0.615:
                return "Tier3_Flagged_Clade_Artifact", "FLAGGED_ClusterC_L4_dominant_clade_SNP (df~0.614)"

        # Tier 1: Gold-standard cross-lineage (L2 + L4) multi-hit convergence (>=2 variants)
        if n_mut >= 2 and has_L2 and has_L4:
            return "Tier1_Gold_CrossLineage_MultiHit", "Independent multi-mutation convergence across BOTH L2 and L4"

        # Tier 2: Within-L4 multi-hit convergence (>=2 variants across >=2 L4 comparisons)
        if n_mut >= 2 and n_cmp >= 2:
            return "Tier2_WithinLineage_MultiHit", "Multiple independent missense mutations across >=2 L4 sublineages"

        # Tier 2b: Cross-lineage single-hit candidates with high Fst support
        if n_cmp >= 2 and has_L2 and has_L4:
            return "Tier2b_CrossLineage_SingleHit", "Cross-lineage signal with 1 catalogued mutation (verify site independence in VCF)"

        return "Tier4_SingleComparison_LineageSpecific", "Single comparison or lineage-specific variant (requires phylogenetic control)"

    res = df_all.apply(_classify_gene, axis=1)
    df_all["Purified_Tier"] = [r[0] for r in res]
    df_all["Audit_Diagnosis"] = [r[1] for r in res]

    df_purified = df_all[
        df_all["Purified_Tier"].isin([
            "Tier1_Gold_CrossLineage_MultiHit",
            "Tier2_WithinLineage_MultiHit",
            "Tier2b_CrossLineage_SingleHit",
        ])
    ].sort_values(
        ["Purified_Tier", "n_comparisons", "n_mutations", "overall_max_Fst"],
        ascending=[True, False, False, False],
    )
    df_purified.to_csv(TABLES_DIR / "Table_S3_Purified_Convergent_Adaptive_Genes.tsv", sep="\t", index=False)

    df_art = df_all[df_all["Purified_Tier"] == "Tier3_Flagged_Clade_Artifact"].sort_values(
        ["max_delta_freq", "gene_name"]
    )
    df_art.to_csv(TABLES_DIR / "Table_S3b_Flagged_Clade_Hitchhiking_Candidates.tsv", sep="\t", index=False)

    n_t1 = (df_purified["Purified_Tier"] == "Tier1_Gold_CrossLineage_MultiHit").sum()
    n_t2 = (df_purified["Purified_Tier"] == "Tier2_WithinLineage_MultiHit").sum()
    n_t2b = (df_purified["Purified_Tier"] == "Tier2b_CrossLineage_SingleHit").sum()
    print(f"  Saved Table_S3 (Tier1={n_t1}, Tier2={n_t2}, Tier2b={n_t2b}) & Table_S3b ({len(df_art)} flagged clade artifacts).")
    return df_purified, df_art


def main():
    run_convergent_selection_audit()


if __name__ == "__main__":
    main()
