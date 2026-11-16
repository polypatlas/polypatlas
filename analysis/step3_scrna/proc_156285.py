#!/usr/bin/env python3
"""Process GSE156285: 4 samples 10x h5 (IT=inferior turbinate, P=polyp; CRS disease atlas, disease-only)."""
import sys, os, gc, json
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad
import scanpy as sc

W = os.path.expanduser("analysis/step3_scrna")
RAW = f"{W}/raw/GSE156285/files"

samples = [("GSM4728360_IT_7741246983_matrix.h5","D1","IT"),
           ("GSM4728361_P_2009946466_matrix.h5","D2","P"),
           ("GSM4728362_1051762513_3976_IT_matrix.h5","D3","IT"),
           ("GSM4728363_1572874946_3976_P_matrix.h5","D4","P")]
ads = []
for fn, donor, tissue in samples:
    a = sc.read_10x_h5(f"{RAW}/{fn}", gex_only=False)
    a.var_names = a.var["gene_symbols"].astype(str).values if "gene_symbols" in a.var else a.var_names
    a.obs_names = [f"{donor}__{b}" for b in a.obs_names]
    a.obs["donor"] = donor; a.obs["tissue"] = tissue; a.obs["group"] = "CRS"
    a = qc_adata(a)
    ads.append(a); print(donor, tissue, a.shape, flush=True); gc.collect()

a = ad.concat(ads); del ads; gc.collect()
print("concat:", a.shape)
ctrl, fails = positive_controls(a)
print("positive controls:", ctrl)
if not (ctrl["TPI1"] > 0.3 and ctrl["GAPDH"] > 0.5):
    print("ASSERT FAIL: housekeeping missing"); sys.exit(1)
print("ASSERT OK (tissue atlas; EPCAM/PTPRC reported, context = mixed tissue)")
a.raw = a
a = norm_log(a); a = cluster(a, res=1.0)
print(a.obs.leiden.value_counts())
sig = signature_scores(a)
pd.set_option("display.width", 200)
print(sig.round(2).to_string())
sig.to_csv(f"{W}/prep/GSE156285_cluster_scores.csv")
a.obs.to_csv(f"{W}/prep/GSE156285_obs.csv")
print("done stage 1")
