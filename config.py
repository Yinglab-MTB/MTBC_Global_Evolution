# -*- coding: utf-8 -*-
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent

DATA_DIR = Path(os.environ.get("MTBC_DATA_DIR", REPO_ROOT / "data")).resolve()
RESULTS_DIR = Path(os.environ.get("MTBC_RESULTS_DIR", REPO_ROOT / "results")).resolve()

TABLES_DIR = RESULTS_DIR / "tables"
ASSETS_DIR = REPO_ROOT / "assets"

SUB14 = [
    "lineage1.1", "lineage1.2", "lineage2.1", "lineage2.2", "lineage3.1",
    "lineage4.1", "lineage4.2", "lineage4.3", "lineage4.4", "lineage4.5",
    "lineage4.6", "lineage4.7", "lineage4.8", "lineage4.9",
]
GLOBAL_SET = {"lineage2.2", "lineage4.1", "lineage4.3", "lineage4.4", "lineage4.8"}
ROOTED_SET = {"lineage2.1", "lineage4.5"}


def ensure_output_dirs():
    for d in (RESULTS_DIR, TABLES_DIR):
        d.mkdir(parents=True, exist_ok=True)
