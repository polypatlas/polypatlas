#!/usr/bin/env python3
"""Stage 2 (generic): concat per-sample h5ad -> positive controls -> cluster -> signature scores -> obs CSV.
Usage: stage2.py <GSE> <res> [extra_obs_col=val ...]
Per-sample h5ad expected in prep/<GSE>_cells/*.h5ad with obs donor/group already set.
Writes prep/<GSE>_cluster_scores.csv and prep/<GSE>_obs.csv. Finalize separately."""
import sys, os, gc
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
from pipeline import *
import anndata as ad, glob

W = os.path.expanduser("analysis/step3_scrna")
GSE, res = sys.argv[1], float(sys.argv[2])
files = sorted(glob.glob(f"{W}/prep/{GSE}_cells/*.h5ad"))
print("samples:", [os.path.basename(f) for f in files])
ads = [ad.read_h5ad(f) for f in files]
a = ad.concat(ads); del ads; gc.collect()
print("concat:", a.shape)
ctrl, fails = positive_controls(a)
print("positive controls:", ctrl)
if not (ctrl.get("TPI1", 0) > 0.3 and ctrl.get("GAPDH", 0) > 0.5):
    print("ASSERT FAIL: housekeeping missing"); sys.exit(1)
print("ASSERT note: EPCAM/PTPRC context judged per dataset in runlog")
a.raw = a
a = norm_log(a); a = cluster(a, res=res)
print(a.obs.leiden.value_counts())
a.obs.to_csv(f"{W}/prep/{GSE}_obs.csv")  # checkpoint before scoring (scoring can OOM on big sets)
print("obs checkpoint saved", flush=True)
sig = signature_scores(a)
pd.set_option("display.width", 220)
print(sig.round(2).to_string())
sig.to_csv(f"{W}/prep/{GSE}_cluster_scores.csv")
a.obs.to_csv(f"{W}/prep/{GSE}_obs.csv")
print("stage 2 done")
