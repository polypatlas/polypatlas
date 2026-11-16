#!/usr/bin/env python3
"""Process GSE276503: 33 samples (HC_IT x2, CRS IT/MT/NP x31). Upgraded pending->core (nasal CRS).
Donor structure unstated -> treat each sample as the donor unit. Stage 1: per-sample QC -> h5ad."""
import sys, os, gc, tarfile, re
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE276503"
OUT = f"{W}/prep/GSE276503_cells"; os.makedirs(OUT, exist_ok=True)

tar_path = f"{RAW}/GSE276503_RAW.tar"
with tarfile.open(tar_path) as tf:
    allf = [m.name for m in tf.getmembers() if m.isfile()]
mtx_files = sorted([f for f in allf if f.endswith("_matrix.mtx.gz")])

def donor_of(fn):
    # GSM8499604_matrix.mtx.gz -> CRSwNP_NP33 (no infix; from GSM titles)
    m = re.search(r"GSM\d+_([A-Za-z0-9]+)_matrix", fn)
    if m: return m.group(1)
    if "GSM8499604_" in fn: return "NP33"
    raise ValueError(f"cannot parse donor from {fn}")

def group_of(d):
    if d.startswith("HC"): return "control"
    return "CRSwNP"

done = 0
for fn in mtx_files:
    donor = donor_of(fn)
    out = f"{OUT}/{donor}.h5ad"
    if os.path.exists(out):
        print(donor, "cached"); continue
    grp = group_of(donor)
    if donor == "NP33":
        f_mtx, f_bc, f_fe = "GSM8499604_matrix.mtx.gz", "GSM8499604_CRSNP33_barcodes.tsv.gz", "GSM8499604_CRSNP33_genes.tsv.gz"
        want = [f_mtx, f_bc, f_fe]
    else:
        base = fn.replace("_matrix.mtx.gz", "")
        want = [n for n in allf if n.startswith(base)]
        f_mtx = base + "_matrix.mtx.gz"
        f_bc = base + "_barcodes.tsv.gz"
        f_cands = [base + "_features.tsv.gz", base + "_genes.tsv.gz"]
        f_fe = next((c for c in f_cands if c in want), None)
    with tarfile.open(tar_path) as tf:
        tf.extractall(f"{RAW}/tmp276", members=[tf.getmember(n) for n in want])
    T = f"{RAW}/tmp276"
    try:
        a = ad.read_mtx(f"{T}/{f_mtx}").T
    except Exception as e:
        print(donor, "SKIP unreadable matrix:", str(e)[:100], flush=True)
        open(f"{W}/prep/{donor}_SKIP.txt","w").write(f"unreadable matrix: {e}")
        for n in want:
            try: os.remove(f"{T}/{n}")
            except OSError: pass
        continue
    fdf = pd.read_csv(f"{T}/{f_fe}", sep="\t", header=None, compression="infer")
    a.var_names = fdf.iloc[:,1].astype(str).values if fdf.shape[1] > 1 else fdf.iloc[:,0].astype(str).values
    a.obs_names = [f"{donor}__{b}" for b in pd.read_csv(f"{T}/{f_bc}", header=None, compression="infer")[0].astype(str).values]
    a.obs["donor"] = donor; a.obs["group"] = grp
    a = qc_adata(a)
    a.write_h5ad(out, compression="gzip")
    print(donor, grp, a.shape, "saved", flush=True)
    del a; gc.collect()
    for n in want:
        try: os.remove(f"{T}/{n}")
        except OSError: pass
    done += 1
print("stage 1 done", done)
