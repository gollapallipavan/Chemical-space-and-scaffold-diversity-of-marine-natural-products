"""Step 2: collapse duplicate 3D-conformer records to one row per COMPOUND_ID.

Keeps the FIRST record per COMPOUND_ID (as stated in Manuscript Section 2.2).
Usage: python 02_deduplicate.py cmnpd_descriptors.csv Supplementary_Data_S1_descriptors_deduplicated.csv
"""
import sys
import pandas as pd

src, dst = sys.argv[1], sys.argv[2]
df = pd.read_csv(src, keep_default_na=False, na_values=[])
dedup = df.drop_duplicates(subset="compound_id", keep="first")
dedup.to_csv(dst, index=False)
print(f"records={len(df):,} unique_compounds={len(dedup):,} "
      f"duplicate_records={len(df)-len(dedup):,} max_records_per_compound={df.compound_id.value_counts().max()}")
