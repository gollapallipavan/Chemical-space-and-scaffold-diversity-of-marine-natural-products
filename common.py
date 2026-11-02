"""Shared paths and helpers for the CMNPD pipeline."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
RES = ROOT / "results"
TABLES = RES / "tables"
FIGS = RES / "figures"
INTER = RES / "intermediate"
for d in (TABLES, FIGS, INTER):
    d.mkdir(parents=True, exist_ok=True)

SEED = 42
S1 = DATA / "Supplementary_Data_S1_descriptors_deduplicated.csv"
S4 = DATA / "Supplementary_Data_S4_full_MCPS_scores.csv"


def read_s1():
    # keep_default_na=False: acyclic compounds have an empty Murcko scaffold string
    # (kept as '' here and counted as one scaffold group, as in the original S4 table).
    return pd.read_csv(S1, keep_default_na=False, na_values=[])
