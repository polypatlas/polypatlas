#!/usr/bin/env python3
"""GSE276503 stage 2 with batched concat (memory discipline: 29 samples)."""
import sys, os, gc, glob
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
files = sorted(glob.glob(f"{W}/prep/GSE276503_cells/*.h5ad"))
print("samples:", [os.path.basename(f) for f in files], flush=True)

# batch concat in groups of 8
parts = []
for i in range(0, len(files), 8):
    chunk = files[i:i+8]
    ads = [ad.read_h5ad(f) for f in chunk]
    c = ad.concat(ads); del ads; gc.collect()
    p = f"{W}/prep/GSE276503_part{i//8}.h5ad"
    c.write_h5ad(p, compression="gzip"); parts.append(p)
    print(f"part {i//8}: {c.shape} saved", flush=True)
    del c; gc.collect()

ads = [ad.read_h5ad(p) for p in parts]
a = ad.concat(ads); del ads; gc.collect()
print("concat:", a.shape, flush=True)
print("positive controls:", positive_controls(a), flush=True)
print("ASSERT note: EPCAM/PTPRC context judged per dataset in runlog", flush=True)
a.raw = a
a = norm_log(a); a = cluster(a, res=1.0)
print(a.obs.leiden.value_counts())
a.obs.to_csv(f"{W}/prep/GSE276503_obs.csv")
print("obs checkpoint saved", flush=True)
sig = signature_scores(a)
pd.set_option("display.width", 220)
print(sig.round(2).to_string())
sig.to_csv(f"{W}/prep/GSE276503_cluster_scores.csv")
a.obs.to_csv(f"{W}/prep/GSE276503_obs.csv")
print("stage 2 done")
