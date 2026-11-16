#!/usr/bin/env python3
"""Process GSE324445: sorted ILC2, NP (6 donors) vs PB (4 donors). Stage 1: per-sample QC -> h5ad."""
import sys, os, gc, json
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE324445/files"
OUT = f"{W}/prep/GSE324445_cells"; os.makedirs(OUT, exist_ok=True)

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
for fn, donor, grp in samples:
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    df = pd.read_csv(f"{RAW}/{fn}", comment="#", index_col=0)
    genes = df.columns.astype(str).values
    X = sp.csr_matrix(df.to_numpy(dtype=np.float32)); del df; gc.collect()
    a = ad.AnnData(X=X)
    a.var_names = genes
    a.obs_names = [f"{donor}__{i}" for i in range(a.n_obs)]
    a.obs["donor"] = donor; a.obs["group"] = grp
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, a.shape, "saved", flush=True)
    del a, X; gc.collect()
print("stage 1 done")
