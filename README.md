# CMNPD cheminformatic profiling — analysis code and data

Code and supplementary data for the manuscript *"Comprehensive Cheminformatic Profiling of the Comprehensive Marine Natural Products Database (CMNPD v1.0): Chemical Space, QSAR Modeling, and a Scaffold-Rarity-Weighted Priority Score for Drug Discovery"* (revised version).

> Replace the placeholders (authors, manuscript DOI/citation, licence holder) before publishing.

## Repository layout

```
data/                         inputs shipped with the manuscript (deduplicated descriptors, MCPS scores, S2, S3)
scripts/
  common.py                   paths, seed (42), shared loader
  01_extract_descriptors.py   RDKit descriptors per SDF record            (needs raw CMNPD SDF)
  02_deduplicate.py           47,451 records -> 25,224 unique COMPOUND_IDs (needs output of 01)
  03_scaffold_mcps.py         Murcko scaffold statistics + MCPS + checks against S2/S3/S4
  04_extended_analysis.py     Sections 3.9-3.13: rule filters, PCA, Chao1, scaffold families, RF triage, cascade
  05_make_figures.py          Figures 6, 7, 8
results/
  tables/                     Supplementary Data S5, S6; Tables S6-S10 (+ regenerated S2/S3/S4)
  figures/                    Figures 6-8
  intermediate/               arrays/tables passed from 04 to 05
  log_*.txt                   console output of the last run
run_all.sh                    runs steps 03-05
```

## Quick start (steps 3-5, no SDF needed)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
./run_all.sh
```

## Full pipeline from the raw SDF (steps 1-2)

Download `CMNPD_1_0_3d_sdf.gz` from <https://www.cmnpd.org> (check the database terms before redistributing it; it is git-ignored here).

```bash
python scripts/01_extract_descriptors.py CMNPD_1_0_3d_sdf.gz cmnpd_descriptors.csv
python scripts/02_deduplicate.py cmnpd_descriptors.csv data/Supplementary_Data_S1_descriptors_deduplicated.csv
./run_all.sh
```

## Reproducibility checks (what was verified)

| Check | Result |
|---|---|
| MCPS regenerated from S1 vs shipped S4 | max abs. difference 2.2e-16 |
| Scaffold frequencies vs S4 | identical |
| Top-50 MCPS IDs vs Table S3 | identical |
| Scaffold counts (7,791 ring scaffolds, 4,336 singletons, 679/3,211 for 50%/80% coverage) | match manuscript |
| QED singleton vs recurring (n = 4,336 vs 20,888) | 0.494 vs 0.473, p = 7.2e-4 (matches) |
| Kruskal-Wallis on MCPS by signature | H = 886.6, p = 7.0e-192 (matches) |
| Tables S5, S6, S6a/b, S7, S8, S9, S10b | numerically identical to the delivered files |
| RF triage classifier | AUC-ROC 0.941, AUC-PR 0.237 (matches, with pinned versions) |

**Version pinning matters.** `scikit-learn==1.8.0` and `rdkit==2026.3.5` reproduce the manuscript numbers exactly. With scikit-learn 1.9.1 the same code gave AUC-PR 0.254 (AUC-ROC 0.942). Over ten random seeds the AUC-PR is 0.223 +/- 0.061, so the difference is within seed noise, but the reported number is only exactly reproducible with the pinned versions.

## Not included in this repository

The following manuscript components were produced by earlier scripts that were not available when this repository was assembled, and were **not** re-created here rather than risk shipping unverified code:

* Record-level analyses of Sections 3.2-3.6 (PCA of ECFP4 bits, K-Means k = 8, drug-likeness / halogenation / novelty Random Forests, Isolation Forest, Williams plot) and Figures 1-3.
* The plotting code for Figures 4 and 5 (their underlying numbers *are* reproduced by `03_scaffold_mcps.py`).

## Citation

[Add manuscript citation / DOI here]. Please also cite the CMNPD paper: Qin JH *et al.* (2021) *Nucleic Acids Res.* 49(D1), D509-D515.
