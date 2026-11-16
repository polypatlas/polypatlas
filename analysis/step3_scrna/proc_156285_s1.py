#!/usr/bin/env python3
"""Process GSE156285: 4 samples 10x h5 (IT/P; CRS atlas). Stage 1: per-sample QC -> h5ad."""
import sys, os, gc
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad
import scanpy as sc

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE156285/files"
OUT = f"{W}/prep/GSE156285_cells"; os.makedirs(OUT, exist_ok=True)

samples = [("GSM4728360_IT_7741246983_matrix.h5","D1","IT"),
           ("GSM4728361_P_2009946466_matrix.h5","D2","P"),
           ("GSM4728362_1051762513_3976_IT_matrix.h5","D3","IT"),
           ("GSM4728363_1572874946_3976_P_matrix.h5","D4","P")]
for fn, donor, tissue in samples:
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    a = sc.read_10x_h5(f"{RAW}/{fn}", gex_only=False)
    a.var_names = a.var["gene_symbols"].astype(str).values if "gene_symbols" in a.var else a.var_names
    a.obs_names = [f"{donor}__{b}" for b in a.obs_names]
    a.obs["donor"] = donor; a.obs["tissue"] = tissue; a.obs["group"] = "CRS"
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, tissue, a.shape, "saved", flush=True)
    del a; gc.collect()
print("stage 1 done")
