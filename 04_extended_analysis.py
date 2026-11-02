"""Step 4: all NEW analyses of Manuscript Sections 3.9-3.13 (Tables 6-10, Supp. S5-S10).

Run after 03. Reads data/S1 and data/S4; writes results/tables/* and results/intermediate/*
(the intermediates are consumed by 05_make_figures.py).
"""
import json
import numpy as np
import pandas as pd
from scipy.stats import entropy, fisher_exact, chi2_contingency
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import average_precision_score, confusion_matrix, precision_recall_curve, roc_auc_score, roc_curve
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem
from common import *

RDLogger.DisableLog("rdApp.*")

df = read_s1()
mc = pd.read_csv(S4, keep_default_na=False, na_values=[])[["compound_id", "scaffold_freq", "MCPS"]]
df = df.merge(mc, on="compound_id", how="left", validate="1:1")
assert df.MCPS.notna().all()

# ---- 2.8 rule-based filters ------------------------------------------------
df["ro5_violations"] = (df.mol_wt > 500).astype(int) + (df.logp > 5) + (df.hbd > 5) + (df.hba > 10)
df["ro5_pass"] = df.ro5_violations <= 1
df["veber_pass"] = (df.rotb <= 10) & (df.tpsa <= 140)
df["ghose_pass"] = df.mol_wt.between(160, 480) & df.logp.between(-0.4, 5.6) & df.mol_mr.between(40, 130) & df.n_heavy_atoms.between(20, 70)
df["egan_pass"] = (df.tpsa <= 131.6) & (df.logp <= 5.88)
df["n_rules_passed"] = df[["ro5_pass", "veber_pass", "ghose_pass", "egan_pass"]].sum(axis=1)


def size_class(r):
    if r.mol_wt <= 300 and r.logp <= 3: return "Fragment-like"
    if r.mol_wt <= 350 and -1 <= r.logp <= 3 and r.rotb <= 7: return "Lead-like"
    return "Drug-like (Ro5)" if r.ro5_pass else "Beyond-Ro5 (bRo5)"


df["size_class"] = df.apply(size_class, axis=1)
df["halogenated"] = df.n_halogens > 0
df["sulfurous"] = df.n_sulfur > 0
df["stereo_density"] = df.n_stereocenters / df.n_heavy_atoms

pd.DataFrame([{"Rule/filter": n, "n_pass": int(df[c].sum()), "pct_pass": round(100 * df[c].mean(), 2)}
              for n, c in [("Lipinski Ro5 (\u22641 violation)", "ro5_pass"), ("Veber rule", "veber_pass"),
                           ("Ghose filter", "ghose_pass"), ("Egan (BOILED-Egg)", "egan_pass")]]
             ).to_csv(TABLES / "Supplementary_Table_S6a_rule_compliance.csv", index=False)
sz = df.size_class.value_counts().rename_axis("size_class").reset_index(name="n_compounds")
sz["pct"] = (100 * sz.n_compounds / len(df)).round(2)
sz.to_csv(TABLES / "Supplementary_Table_S6b_size_classification.csv", index=False)
s5 = ["compound_id", "mol_wt", "logp", "tpsa", "hbd", "hba", "rotb", "n_rings", "n_arom_rings", "frac_csp3", "mol_mr",
      "n_stereocenters", "n_heavy_atoms", "stereo_density", "n_halogens", "n_sulfur", "qed", "murcko_scaffold_smiles",
      "scaffold_freq", "MCPS", "ro5_violations", "ro5_pass", "veber_pass", "ghose_pass", "egan_pass", "n_rules_passed", "size_class"]
df[s5].to_csv(TABLES / "Supplementary_Data_S5_rule_compliance_stereo_density.csv", index=False)

# ---- 2.9 correlation + descriptor PCA --------------------------------------
desc = ["mol_wt", "logp", "tpsa", "hbd", "hba", "rotb", "n_rings", "n_arom_rings", "frac_csp3", "mol_mr",
        "n_stereocenters", "n_heavy_atoms", "qed"]
df[desc].corr().to_csv(INTER / "descriptor_correlation.csv")
pca = PCA(n_components=5, random_state=SEED).fit(StandardScaler().fit_transform(df[desc]))
pd.DataFrame(pca.components_.T, index=desc, columns=[f"PC{i+1}" for i in range(5)]).to_csv(INTER / "pca_loadings.csv")
print("descriptor PCA explained variance:", np.round(pca.explained_variance_ratio_[:2], 4))

# ---- 2.10 scaffold richness / diversity ------------------------------------
ring = df[df.murcko_scaffold_smiles != ""]
counts = ring.murcko_scaffold_smiles.value_counts()
S_obs, N = len(counts), len(ring)
f1, f2 = int((counts == 1).sum()), int((counts == 2).sum())
chao1 = S_obs + f1 * (f1 - 1) / (2 * (f2 + 1))
se = np.sqrt(f2 * (0.5 * (f1 / f2) ** 2 + (f1 / f2) ** 3 + 0.25 * (f1 / f2) ** 4))
p = counts.values / N
H = entropy(p)
cum = np.cumsum(counts.values) / N * 100
stats = dict(S_obs=S_obs, N=N, f1=f1, f2=f2, chao1=chao1, chao1_se=se, shannon_nats=H, shannon_bits=H / np.log(2),
             simpson=1 - (p ** 2).sum(), inv_simpson=1 / (p ** 2).sum(), pielou=H / np.log(S_obs),
             n50=int(np.searchsorted(cum, 50) + 1), n80=int(np.searchsorted(cum, 80) + 1))
pd.DataFrame([stats]).to_csv(TABLES / "Supplementary_Table_S7_scaffold_diversity_indices.csv", index=False)
json.dump({k: float(v) for k, v in stats.items()}, open(INTER / "scaffold_diversity_stats.json", "w"), indent=2)

rng = np.random.default_rng(SEED)
arr = ring.murcko_scaffold_smiles.values
sizes = np.unique(np.round(np.logspace(np.log10(50), np.log10(N), 40)).astype(int))
curve = np.zeros((15, len(sizes)))
for r in range(15):
    sh, seen, ptr = arr[rng.permutation(N)], set(), 0
    for j, n in enumerate(sizes):
        seen.update(sh[ptr:n]); ptr = n; curve[r, j] = len(seen)
pd.DataFrame({"sample_size": sizes, "mean_unique_scaffolds": curve.mean(0), "std": curve.std(0)}
             ).to_csv(INTER / "rarefaction_curve.csv", index=False)

# ---- scaffold fingerprints, PCA, K-Means families --------------------------
agg = (ring.groupby("murcko_scaffold_smiles")
       .agg(n_compounds=("compound_id", "count"), mean_qed=("qed", "mean"),
            pct_halogenated=("halogenated", lambda x: 100 * x.mean()), pct_sulfur=("sulfurous", lambda x: 100 * x.mean()),
            mean_stereocenters=("n_stereocenters", "mean"), mean_mw=("mol_wt", "mean"), mean_mcps=("MCPS", "mean"))
       .reset_index().rename(columns={"murcko_scaffold_smiles": "scaffold_smiles"})
       .sort_values("n_compounds", ascending=False).reset_index(drop=True))  # default sort: row order (hence K-Means init) matches the manuscript run
fps = np.zeros((len(agg), 512), dtype=np.int8)
for i, smi in enumerate(agg.scaffold_smiles):
    bv = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(smi), 2, nBits=512)
    DataStructs.ConvertToNumpyArray(bv, fps[i])
xy = PCA(n_components=2, random_state=SEED).fit_transform(fps)
agg["pc1"], agg["pc2"] = xy[:, 0], xy[:, 1]
agg["scaffold_family"] = KMeans(n_clusters=10, n_init=10, random_state=SEED).fit_predict(fps)
agg.to_csv(TABLES / "Supplementary_Data_S6_scaffold_family_assignments.csv", index=False)
agg.to_csv(INTER / "scaffold_aggregates_pca.csv", index=False)
fam = (agg.groupby("scaffold_family")
       .apply(lambda g: pd.Series({"n_scaffolds": len(g), "n_compounds": g.n_compounds.sum(),
                                   "wmean_qed": np.average(g.mean_qed, weights=g.n_compounds),
                                   "wmean_halogen_pct": np.average(g.pct_halogenated, weights=g.n_compounds),
                                   "wmean_mcps": np.average(g.mean_mcps, weights=g.n_compounds)}), include_groups=False)
       .reset_index().sort_values("n_compounds", ascending=False))
fam.to_csv(TABLES / "Supplementary_Table_S8_scaffold_family_summary.csv", index=False)

# ---- 2.11 halogen x sulfur, triage RF, cascade -----------------------------
ct = pd.crosstab(df.halogenated, df.sulfurous)
print(f"Fisher OR={fisher_exact(ct)[0]:.3f} p={fisher_exact(ct)[1]:.4f}; chi2={chi2_contingency(ct)[0]:.2f} p={chi2_contingency(ct)[1]:.4f}")

t99, t90 = df.MCPS.quantile(0.99), df.MCPS.quantile(0.90)
df["top1pct"] = (df.MCPS >= t99).astype(int)
feats = ["mol_wt", "logp", "tpsa", "hbd", "hba", "rotb", "n_rings", "n_arom_rings", "frac_csp3", "mol_mr",
         "n_stereocenters", "n_heavy_atoms"]  # QED, scaffold frequency, halogen/sulfur deliberately excluded
Xtr, Xte, ytr, yte = train_test_split(df[feats].values, df.top1pct.values, test_size=0.2, random_state=SEED, stratify=df.top1pct.values)
rf = RandomForestClassifier(n_estimators=300, max_depth=8, class_weight="balanced", random_state=SEED, n_jobs=-1).fit(Xtr, ytr)
proba = rf.predict_proba(Xte)[:, 1]
auc, ap = roc_auc_score(yte, proba), average_precision_score(yte, proba)
cm = confusion_matrix(yte, proba >= 0.5)
fpr, tpr, _ = roc_curve(yte, proba)
prec, rec, _ = precision_recall_curve(yte, proba)
np.savez(INTER / "rf_triage_results.npz", fpr=fpr, tpr=tpr, prec=prec, rec=rec, cm=cm, auc=auc, ap=ap, yte=yte, proba=proba)
imp = pd.Series(rf.feature_importances_, index=feats).sort_values(ascending=False)
imp.to_csv(TABLES / "Supplementary_Table_S10b_RF_triage_feature_importance.csv")
pd.DataFrame([{"model": "RandomForestClassifier(n_estimators=300, max_depth=8, class_weight=balanced)",
               "features": ", ".join(feats), "target": "MCPS \u2265 99th percentile (top-1% tier)",
               "train_test_split": "80/20 stratified", "n_test": len(yte), "n_test_positive": int(yte.sum()),
               "test_AUC_ROC": auc, "test_AUC_PR": ap, "baseline_prevalence": yte.mean(),
               "precision_at_0.5": cm[1, 1] / (cm[1, 1] + cm[0, 1]), "recall_at_0.5": cm[1, 1] / (cm[1, 1] + cm[1, 0])}]
             ).to_csv(TABLES / "Supplementary_Table_S10_RF_triage_performance.csv", index=False)
print(f"RF triage AUC-ROC={auc:.3f} AUC-PR={ap:.3f}; confusion={cm.tolist()}")

casc = [len(df), int(df.ro5_pass.sum()), int((df.ro5_pass & df.veber_pass).sum()),
        int((df.ro5_pass & df.veber_pass & (df.MCPS >= t90)).sum()), int((df.ro5_pass & df.veber_pass & (df.MCPS >= t99)).sum())]
print("cascade:", casc)
final = df[df.ro5_pass & df.veber_pass & (df.MCPS >= t99)].sort_values("MCPS", ascending=False)
final[["compound_id", "mol_wt", "logp", "tpsa", "qed", "n_stereocenters", "n_halogens", "n_sulfur",
       "murcko_scaffold_smiles", "scaffold_freq", "MCPS", "ro5_violations"]
      ].to_csv(TABLES / "Supplementary_Table_S9_final_shortlist_252.csv", index=False)
json.dump({"cascade": casc, "t90": float(t90), "t99": float(t99)}, open(INTER / "cascade.json", "w"))
df.to_csv(INTER / "master_table.csv", index=False)
