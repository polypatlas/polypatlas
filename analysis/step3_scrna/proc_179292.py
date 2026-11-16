#!/usr/bin/env python3
"""Process GSE179292: 4 donors, single T-helper cells from nasal polyps (disease-only)."""
import sys, os, gc
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE179292/files"
PY = f"{W}/venv/bin/python"

donors = [("GSM5413781","NP_1","GSM5413781_NP_1_GE"),
          ("GSM5413782","NP_2","GSM5413782_NP_2"),
          ("GSM5413783","NP_3","GSM5413783_NP_3"),
          ("GSM5413784","NP_4","GSM5413784_NP_4")]
ads = []
for gsm, donor, prefix in donors:
    a = ad.read_mtx(f"{RAW}/{prefix}_matrix.mtx.gz").T
    feats = pd.read_csv(f"{RAW}/{prefix}_features.tsv.gz", sep="\t", header=None)
    a.var_names = feats.iloc[:,1].astype(str).values if feats.shape[1] > 1 else feats.iloc[:,0].astype(str).values
    bc = pd.read_csv(f"{RAW}/{prefix}_barcodes.tsv.gz", header=None)[0].astype(str).values
    a.obs_names = [f"{donor}_{b}" for b in bc]
    a.obs["donor"] = donor; a.obs["group"] = "CRSwNP"
    a = qc_adata(a)
    ads.append(a)
    print(donor, a.shape, flush=True)

a = ad.concat(ads, label="donor_src", index_unique="-")
del ads; gc.collect()
print("concat:", a.shape)

ctrl, fails = positive_controls(a)
print("positive controls:", ctrl)
# Context: sorted T-helper prep -> EPCAM absent is biological, not a broken object.
# Object is sound if housekeeping robust + PTPRC dominant.
if not (ctrl["TPI1"] > 0.3 and ctrl["GAPDH"] > 0.5 and ctrl["PTPRC"] > 0.5):
    print("ASSERT FAIL: housekeeping/immune signal missing"); sys.exit(1)
print("ASSERT OK (sorted T-cell prep; EPCAM absent as expected; note written to runlog)")

a.raw = a
a = norm_log(a)
a = cluster(a)
print(a.obs.leiden.value_counts())
sig = signature_scores(a)
pd.set_option("display.width", 200)
print(sig.round(2).to_string())
sig.to_csv(f"{W}/prep/GSE179292_cluster_scores.csv")

# mapping written after inspecting scores
mapping = {}  # fill below
print("WRITE MAPPING MANUALLY AFTER INSPECTION")
