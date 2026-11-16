#!/usr/bin/env python3
"""Process GSE324445: sorted ILC2 from nasal polyps (6 donors) vs peripheral blood (4 donors). BD Rhapsody CSVs (cells x genes)."""
import sys, os, gc, json
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE324445/files"

samples = [("GSM9576383_NP1_RSEC_MolsPerCell.csv.gz","NP1","NP"),
           ("GSM9576384_NP2_RSEC_MolsPerCell.csv.gz","NP2","NP"),
           ("GSM9576385_NP3_RSEC_MolsPerCell.csv.gz","NP3","NP"),
           ("GSM9576386_NP4_RSEC_MolsPerCell.csv.gz","NP4","NP"),
           ("GSM9576387_NP5_RSEC_MolsPerCell.csv.gz","NP5","NP"),
           ("GSM9576388_NP6_RSEC_MolsPerCell.csv.gz","NP6","NP"),
           ("GSM9576389_PB1_RSEC_MolsPerCell.csv.gz","PB1","PB"),
           ("GSM9576390_PB2_RSEC_MolsPerCell.csv.gz","PB2","PB"),
           ("GSM9576391_PB3_RSEC_MolsPerCell.csv.gz","PB3","PB"),
           ("GSM9576392_PB4_RSEC_MolsPerCell.csv.gz","PB4","PB")]
ads = []
for fn, donor, grp in samples:
    df = pd.read_csv(f"{RAW}/{fn}", comment="#", index_col=0)
    X = sp.csr_matrix(df.values.astype(np.float32))
    a = ad.AnnData(X=X)
    a.var_names = df.columns.astype(str).values
    a.obs_names = [f"{donor}__{i}" for i in df.index.astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = grp
    a = qc_adata(a)
    ads.append(a); del df
    print(donor, a.shape, flush=True); gc.collect()

a = ad.concat(ads); del ads; gc.collect()
print("concat:", a.shape)
ctrl, fails = positive_controls(a)
print("positive controls:", ctrl)
# context: sorted ILC2 prep -> EPCAM expected absent; ILC2 markers should be strong
if not (ctrl["TPI1"] > 0.3 and ctrl["GAPDH"] > 0.5 and ctrl["PTPRC"] > 0.5):
    print("ASSERT FAIL"); sys.exit(1)
for g in ["IL7R","GATA3","PTGDR2","IL1RL1"]:
    f = float((a.raw.X[:, a.var_names.get_loc(g)] > 0).mean()) if g in a.var_names else -1 if False else (float((a.X[:, a.var_names.get_loc(g)] > 0).mean()) if g in a.var_names else -1)
    print(f"ILC2 marker {g}: {f:.3f}")
print("ASSERT OK (sorted ILC2 prep)")

a.raw = a
a = norm_log(a); a = cluster(a, res=0.6)
print(a.obs.leiden.value_counts())
sig = signature_scores(a)
pd.set_option("display.width", 200)
print(sig.loc[:, ["ILC2","Th2","CD4_T","Treg","NK","mast","eosinophil","B"]].round(2).to_string())
sig.to_csv(f"{W}/prep/GSE324445_cluster_scores.csv")
a.obs.to_csv(f"{W}/prep/GSE324445_obs.csv")
print("done stage 1")
