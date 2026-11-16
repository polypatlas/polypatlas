#!/usr/bin/env python3
"""Process GSE202100: 5 nasal polyp donors (disease-only, PMID 35608904, tuft cell study).
Stage 1: per-sample QC -> h5ad."""
import sys, os, gc, tarfile
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE202100"
OUT = f"{W}/prep/GSE202100_cells"; os.makedirs(OUT, exist_ok=True)
T = f"{RAW}/tmp"; os.makedirs(T, exist_ok=True)

donors = [("GSM6094585_Polyp1","Polyp1"),("GSM6094586_Polyp2","Polyp2"),("GSM6094587_Polyp3","Polyp3"),
          ("GSM6094588_Polyp4","Polyp4"),("GSM6094589_Polyp5","Polyp5")]
tar_path = f"{RAW}/GSE202100_RAW.tar"
with tarfile.open(tar_path) as tf:
    allf = [m.name for m in tf.getmembers() if m.isfile()]

for prefix, donor in donors:
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    want = [n for n in allf if n.startswith(prefix)]
    with tarfile.open(tar_path) as tf:
        tf.extractall(T, members=[tf.getmember(n) for n in want])
    a = ad.read_mtx(f"{T}/{prefix}_matrix.mtx.gz").T
    fdf = pd.read_csv(f"{T}/{prefix}_features.tsv.gz", sep="\t", header=None, compression="infer")
    a.var_names = fdf.iloc[:,1].astype(str).values if fdf.shape[1] > 1 else fdf.iloc[:,0].astype(str).values
    a.obs_names = [f"{donor}__{b}" for b in pd.read_csv(f"{T}/{prefix}_barcodes.tsv.gz", header=None, compression="infer")[0].astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = "CRSwNP"
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, a.shape, "saved", flush=True)
    del a; gc.collect()
    for n in want:
        try: os.remove(f"{T}/{n}")
        except OSError: pass
print("stage 1 done")
