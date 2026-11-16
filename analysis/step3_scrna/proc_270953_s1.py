#!/usr/bin/env python3
"""Process GSE270953 scRNA arm: B lineage, nasal polyp (5 donors) vs tonsil (5 donors). Stage 1: per-sample QC -> h5ad."""
import sys, os, gc
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE270953/files"
OUT = f"{W}/prep/GSE270953_cells"; os.makedirs(OUT, exist_ok=True)

np_donors = ["NP21009","NP21010","NP21040","NP22069","NP22077"]
tn_donors = ["TN19056","TN21006","TN21055","TN22071","TN22072"]
for donor in np_donors + tn_donors:
    grp = "NP" if donor.startswith("NP") else "tonsil"
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    pfx = f"GSE270953_{donor}"
    a = ad.read_mtx(f"{RAW}/{pfx}_matrix.mtx.gz").T
    f = pd.read_csv(f"{RAW}/{pfx}_features.tsv.gz", sep="\t", header=None)
    a.var_names = f.iloc[:,1].astype(str).values if f.shape[1] > 1 else f.iloc[:,0].astype(str).values
    a.obs_names = [f"{donor}__{b}" for b in pd.read_csv(f"{RAW}/{pfx}_barcodes.tsv.gz", header=None)[0].astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = grp
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, grp, a.shape, "saved", flush=True)
    del a; gc.collect()
print("stage 1 done")
