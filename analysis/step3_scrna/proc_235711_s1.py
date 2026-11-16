#!/usr/bin/env python3
"""Process GSE235711 (Liao): 20 donors, per-donor filtered.tar.gz (each a tar of barcodes/features/matrix).
Groups: Control / CRSsNP_Eth / CRSwNP_Eth / CRSwNP_NP. Stage 1: per-sample QC -> h5ad."""
import sys, os, gc, tarfile, re, shutil
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE235711"
OUT = f"{W}/prep/GSE235711_cells"; os.makedirs(OUT, exist_ok=True)
TMP = f"{RAW}/tmp_inner"; os.makedirs(TMP, exist_ok=True)

def group_of(donor):
    if donor.startswith("Control"): return "control"
    if donor.startswith("CRSsNP"): return "CRSsNP"
    if donor.startswith("CRSwNP"): return "CRSwNP"
    return "unknown"

tar_path = f"{RAW}/GSE235711_RAW.tar"
with tarfile.open(tar_path) as tf:
    members = [m.name for m in tf.getmembers() if m.isfile()]

for mem in sorted(members):
    donor = re.sub(r"^GSM\d+_", "", mem).replace("_filtered.tar.gz", "")
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    grp = group_of(donor)
    # extract inner tar
    for f_ in os.listdir(TMP):
        p = os.path.join(TMP, f_)
        os.remove(p) if os.path.isfile(p) else shutil.rmtree(p)
    with tarfile.open(tar_path) as tf:
        inner_data = tf.extractfile(mem).read()
    inner_path = os.path.join(TMP, "inner.tar.gz")
    open(inner_path, "wb").write(inner_data)
    with tarfile.open(inner_path) as itf:
        itf.extractall(TMP)
    a = ad.read_mtx(f"{TMP}/matrix.mtx.gz").T
    fdf = pd.read_csv(f"{TMP}/features.tsv.gz", sep="\t", header=None, compression="infer")
    a.var_names = fdf.iloc[:,1].astype(str).values if fdf.shape[1] > 1 else fdf.iloc[:,0].astype(str).values
    a.obs_names = [f"{donor}__{b}" for b in pd.read_csv(f"{TMP}/barcodes.tsv.gz", header=None, compression="infer")[0].astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = grp
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, grp, a.shape, "saved", flush=True)
    del a; gc.collect()
shutil.rmtree(TMP, ignore_errors=True)
print("stage 1 done")
