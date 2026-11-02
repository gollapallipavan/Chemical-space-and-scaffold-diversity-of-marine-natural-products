"""Step 3: Bemis-Murcko scaffold statistics and the Marine Chemotype Priority Score (MCPS).

MCPS = 0.40*z(QED) + 0.40*z(-log10 scaffold_frequency) + 0.20*z(marine_signature)
marine_signature = 0.6*[halogenated] + 0.4*[sulfur-containing]
z() = standard score with the sample SD (ddof=1) across all 25,224 compounds.

Outputs are checked against the shipped Supplementary Data S4 / Tables S2, S3.
"""
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, kruskal
from common import *

df = read_s1()
freq = df["murcko_scaffold_smiles"].map(df["murcko_scaffold_smiles"].value_counts())
df["scaffold_freq"] = freq  # acyclic compounds (scaffold == '') form one group of 1,098

z = lambda x: (x - x.mean()) / x.std(ddof=1)
sig = 0.6 * (df.n_halogens > 0) + 0.4 * (df.n_sulfur > 0)
df["MCPS"] = 0.4 * z(df.qed) + 0.4 * z(-np.log10(df.scaffold_freq)) + 0.2 * z(sig)

cols = ["compound_id", "mol_wt", "logp", "tpsa", "qed", "n_stereocenters", "n_halogens",
        "n_sulfur", "scaffold_freq", "murcko_scaffold_smiles", "MCPS"]
df[cols].to_csv(TABLES / "Supplementary_Data_S4_full_MCPS_scores.regenerated.csv", index=False)
df.sort_values("MCPS", ascending=False).head(50)[cols].to_csv(
    TABLES / "Supplementary_Table_S3_top50_MCPS_candidates.regenerated.csv", index=False)

# ---- verification against the shipped S4 / S3 -------------------------------
ref = pd.read_csv(S4, keep_default_na=False, na_values=[])
chk = df.merge(ref[["compound_id", "MCPS", "scaffold_freq"]], on="compound_id", suffixes=("", "_ref"))
print("max |MCPS - S4 MCPS|        :", float(np.abs(chk.MCPS - chk.MCPS_ref).max()))
print("scaffold_freq identical      :", bool((chk.scaffold_freq == chk.scaffold_freq_ref).all()))
s3 = pd.read_csv(DATA / "Supplementary_Table_S3_top50_MCPS_candidates.csv")
print("top-50 IDs identical to S3   :", list(df.nlargest(50, "MCPS").compound_id) == list(s3.compound_id))

# ---- scaffold statistics (Manuscript Section 3.7) ---------------------------
ring = df[df.murcko_scaffold_smiles != ""]
vc = ring.murcko_scaffold_smiles.value_counts()
print(f"unique ring scaffolds={len(vc):,}; singletons={(vc==1).sum():,}; ring compounds={len(ring):,}")
cum = np.cumsum(vc.values) / len(ring) * 100
print("scaffolds for 50% / 80%:", int(np.searchsorted(cum, 50) + 1), int(np.searchsorted(cum, 80) + 1))
# "recurring" = every compound whose scaffold group has >1 member (includes the 1,098 acyclic
# compounds), i.e. n = 20,888, matching Manuscript Section 3.7 / Figure 4D.
single = df.scaffold_freq == 1
u, p = mannwhitneyu(df.qed[single], df.qed[~single], alternative="two-sided")
print(f"QED singleton (n={single.sum():,}) vs recurring (n={(~single).sum():,}): "
      f"median {df.qed[single].median():.3f} vs {df.qed[~single].median():.3f}, MWU p={p:.2e}")

hal, sul = df.n_halogens > 0, df.n_sulfur > 0
grp = np.select([hal & sul, hal, sul], [3, 2, 1], 0)
H, pk = kruskal(*[df.MCPS[grp == g] for g in range(4)])
print(f"Kruskal-Wallis H={H:.1f}, p={pk:.1e}; top-1% threshold={df.MCPS.quantile(0.99):.3f}")

top = (ring.groupby("murcko_scaffold_smiles")
       .agg(n_compounds=("compound_id", "count"), mean_QED=("qed", "mean"),
            pct_halogenated=("n_halogens", lambda x: 100 * (x > 0).mean()),
            mean_stereocenters=("n_stereocenters", "mean"), mean_MW=("mol_wt", "mean"))
       .sort_values("n_compounds", ascending=False).head(20).reset_index())
top.insert(2, "pct_of_dataset", (100 * top.n_compounds / len(df)).round(2))
top.to_csv(TABLES / "Supplementary_Table_S2_top_scaffolds.regenerated.csv", index=False)
