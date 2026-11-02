"""Step 1 (original Supplementary Code S1, only the input/output paths made configurable).
Usage: python 01_extract_descriptors.py CMNPD_1_0_3d_sdf.gz cmnpd_descriptors.csv
Record level: 47,451 SDF records -> one descriptor row per record."""
import gzip, time, csv, sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski, QED
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

t0 = time.time()
path = sys.argv[1] if len(sys.argv) > 1 else 'CMNPD_1_0_3d_sdf.gz'
supplier = Chem.ForwardSDMolSupplier(gzip.open(path, 'rb'))

out = open(sys.argv[2] if len(sys.argv) > 2 else 'cmnpd_descriptors.csv', 'w', newline='')
w = csv.writer(out)
w.writerow(['compound_id','mol_wt','logp','tpsa','hbd','hba','rotb','n_rings','n_arom_rings',
            'frac_csp3','mol_mr','n_stereocenters','n_halogens','n_F','n_Cl','n_Br','n_I',
            'n_sulfur','qed','murcko_scaffold_smiles','n_heavy_atoms'])

n_total = 0
n_valid = 0
n_failed = 0
for mol in supplier:
    n_total += 1
    if mol is None:
        n_failed += 1
        continue
    try:
        cid = mol.GetProp('COMPOUND_ID') if mol.HasProp('COMPOUND_ID') else f'mol_{n_total}'
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rotb = Descriptors.NumRotatableBonds(mol)
        nrings = Lipinski.RingCount(mol)
        narom = Lipinski.NumAromaticRings(mol)
        fcsp3 = Lipinski.FractionCSP3(mol)
        mr = Crippen.MolMR(mol)
        stereo = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False))
        halo_counts = {'F':0,'Cl':0,'Br':0,'I':0}
        n_s = 0
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym in halo_counts:
                halo_counts[sym]+=1
            if sym == 'S':
                n_s += 1
        n_halo = sum(halo_counts.values())
        qed = QED.qed(mol) if mol.GetNumHeavyAtoms() <= 200 else None
        try:
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        except Exception:
            scaf = ''
        nha = mol.GetNumHeavyAtoms()
        w.writerow([cid, round(mw,3), round(logp,3), round(tpsa,3), hbd, hba, rotb, nrings, narom,
                    round(fcsp3,3), round(mr,3), stereo, n_halo, halo_counts['F'], halo_counts['Cl'],
                    halo_counts['Br'], halo_counts['I'], n_s, round(qed,4) if qed is not None else '',
                    scaf, nha])
        n_valid += 1
    except Exception as e:
        n_failed += 1
    if n_total % 5000 == 0:
        print(f'{n_total} processed, {n_valid} valid, {time.time()-t0:.1f}s', file=sys.stderr)

out.close()
print(f'DONE. total={n_total} valid={n_valid} failed={n_failed} time={time.time()-t0:.1f}s')
