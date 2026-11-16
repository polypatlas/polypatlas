#!/usr/bin/env python3
"""Process GSE255544: 7 CRSwNP donors, lymphocytes (donors 1-3 mtx triplets; donors 4-7 pooled h5).
Upgraded pending->core: GSM titles confirm CRSwNP. Disease-only -> composition-only. Stage 1."""
import sys, os, gc
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad
import scanpy as sc

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE255544/files"
OUT = f"{W}/prep/GSE255544_cells"; os.makedirs(OUT, exist_ok=True)

triplets = [("GSM8074487_donor1","donor1"),("GSM8074488_donor2","donor2"),("GSM8074489_donor3","donor3")]
for prefix, donor in triplets:
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    a = ad.read_mtx(f"{RAW}/{prefix}matrix.mtx.gz").T
    f = pd.read_csv(f"{RAW}/{prefix}features.tsv.gz", sep="\t", header=None)
    a.var_names = f.iloc[:,1].astype(str).values if f.shape[1] > 1 else f.iloc[:,0].astype(str).values
    a.obs_names = [f"{donor}__{b}" for b in pd.read_csv(f"{RAW}/{prefix}barcodes.tsv.gz", header=None)[0].astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = "CRSwNP"
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, a.shape, "saved", flush=True); del a; gc.collect()

out = f"{OUT}/donors4-7.h5ad"
if not os.path.exists(out):
    a = sc.read_10x_h5(f"{RAW}/GSM8074490_donors4-7filtered_feature_bc_matrix.h5", gex_only=False)
    a.var_names = a.var["gene_symbols"].astype(str).values if "gene_symbols" in a.var else a.var_names
    a.obs_names = [f"donors4-7__{b}" for b in a.obs_names]
    a.obs["donor"] = "donors4-7"; a.obs["group"] = "CRSwNP"
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print("donors4-7", a.shape, "saved", flush=True); del a; gc.collect()
else:
    print("donors4-7 cached")
print("stage 1 done")
