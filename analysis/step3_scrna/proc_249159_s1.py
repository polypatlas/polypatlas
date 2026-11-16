#!/usr/bin/env python3
"""Process GSE249159 human arm only: ECRS polyp (3) vs PBMC (3), ST2hi memory Th2 sorted.
Files: GSM{g}_{i} where i=5..10 (human). Stage 1: per-sample QC -> h5ad."""
import sys, os, gc
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE249159/files"
OUT = f"{W}/prep/GSE249159_cells"; os.makedirs(OUT, exist_ok=True)

samples = [("GSM7927627_5","PBMC_1","PBMC"),("GSM7927628_6","PBMC_2","PBMC"),
           ("GSM7927629_7","PBMC_3","PBMC"),("GSM7927630_8","Polyp_1","Polyp"),
           ("GSM7927631_9","Polyp_2","Polyp"),("GSM7927632_10","Polyp_3","Polyp")]
for prefix, donor, grp in samples:
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    a = ad.read_mtx(f"{RAW}/{prefix}_matrix.mtx.gz").T
    f = pd.read_csv(f"{RAW}/{prefix}_features.tsv.gz", sep="\t", header=None)
    a.var_names = f.iloc[:,1].astype(str).values if f.shape[1] > 1 else f.iloc[:,0].astype(str).values
    a.obs_names = [f"{donor}__{b}" for b in pd.read_csv(f"{RAW}/{prefix}_barcodes.tsv.gz", header=None)[0].astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = grp
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, grp, a.shape, "saved", flush=True)
    del a; gc.collect()
print("stage 1 done")
