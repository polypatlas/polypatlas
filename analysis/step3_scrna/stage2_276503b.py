#!/usr/bin/env python3
"""GSE276503 stage 2b: cluster each part separately (memory discipline), same params; annotate later by harmonized markers."""
import sys, os, gc, glob
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
parts = sorted(glob.glob(f"{W}/prep/GSE276503_part*.h5ad"))
for p in parts:
    tag = os.path.basename(p).replace(".h5ad", "")
    obs_out = f"{W}/prep/{tag}_obs.csv"
    if os.path.exists(obs_out):
        print(tag, "cached"); continue
    a = ad.read_h5ad(p)
    print(tag, a.shape, flush=True)
    print("positive controls:", positive_controls(a), flush=True)
    a.raw = a
    a = norm_log(a); a = cluster(a, res=1.0)
    print(a.obs.leiden.value_counts())
    a.obs.to_csv(obs_out)
    print(tag, "obs saved", flush=True)
    sig = signature_scores(a)
    sig.to_csv(f"{W}/prep/{tag}_cluster_scores.csv")
    print(tag, "scores saved", flush=True)
    del a; gc.collect()
print("stage 2b done")
