#!/usr/bin/env python3
"""GSE276503 finalize: annotate 4 parts with harmonized labels, save part objects,
combined composition + pseudobulk. Design: control n=2 vs CRSwNP n=27 -> composition-only."""
import sys, os, gc, glob, json
sys.path.insert(0, os.path.expanduser("analysis/step3_scrna"))
import pandas as pd, numpy as np, anndata as ad
from pipeline import composition_table, pseudobulk

W = os.path.expanduser("analysis/step3_scrna")
MAP = {
 "0": {"0":"glandular","1":"CD8_T","2":"fibroblast","3":"ambiguous","4":"endothelial","5":"CD4_T",
       "6":"endothelial","7":"fibroblast","8":"glandular","9":"B","10":"basal","11":"basal",
       "12":"goblet","13":"basal","14":"secretory","15":"glandular","16":"macrophage",
       "17":"ambiguous","18":"mast","19":"fibroblast","20":"basal"},
 "1": {"0":"CD8_T","1":"NKT","2":"glandular","3":"ambiguous","4":"endothelial","5":"CD8_T",
       "6":"goblet","7":"fibroblast","8":"CD8_T","9":"fibroblast","10":"endothelial","11":"mast",
       "12":"B","13":"basal","14":"basal","15":"plasma","16":"macrophage","17":"goblet",
       "18":"ciliated","19":"ambiguous"},
 "2": {"0":"CD8_T","1":"basal","2":"endothelial","3":"fibroblast","4":"basal","5":"B",
       "6":"ambiguous","7":"basal","8":"fibroblast","9":"monocyte","10":"NK","11":"goblet",
       "12":"mast","13":"glandular","14":"CD8_T","15":"basal","16":"basal","17":"plasma",
       "18":"basal","19":"plasma","20":"plasma","21":"ILC2","22":"plasma","23":"pDC","24":"ambiguous"},
 "3": {"0":"glandular","1":"CD8_T","2":"CD4_T","3":"fibroblast","4":"fibroblast","5":"ambiguous",
       "6":"CD8_T","7":"monocyte","8":"NK","9":"basal","10":"B","11":"endothelial","12":"secretory",
       "13":"fibroblast","14":"mast","15":"basal","16":"glandular","17":"plasma","18":"plasma",
       "19":"glandular","20":"ciliated","21":"secretory","22":"basal","23":"plasma","24":"plasma",
       "25":"pDC","26":"ambiguous","27":"ionocyte"},
}
def group_of(d): return "control" if d.startswith("HC") else "CRSwNP"

comp_all, pb_all = [], []
total_cells = 0
for pn, mapping in sorted(MAP.items()):
    p = f"{W}/prep/GSE276503_part{pn}.h5ad"
    obs = pd.read_csv(f"{W}/prep/GSE276503_part{pn}_obs.csv", index_col=0)
    a = ad.read_h5ad(p)
    a.obs["leiden"] = a.obs_names.map(dict(zip(obs.index, obs.leiden.astype(str)))).values
    a.obs["celltype"] = a.obs["leiden"].map(mapping).values
    a.obs["donor"] = a.obs_names.str.split("__").str[0].values
    a.obs["group"] = a.obs["donor"].map(group_of).values
    a.obs["part"] = pn
    print(f"part{pn} unmapped:", int(a.obs.celltype.isna().sum()), flush=True)
    a.write_h5ad(f"{W}/objects/GSE276503_part{pn}_annotated.h5ad", compression="gzip")
    print(f"part{pn} annotated object saved", flush=True)
    total_cells += a.n_obs
    ct = composition_table(a); ct["dataset"] = "GSE276503"
    comp_all.append(ct)
    pb = pseudobulk(a); pb_all.append(pb)
    del a; gc.collect()

pd.concat(comp_all).to_csv(f"{W}/prep/GSE276503_composition.csv", index=False)
pd.concat(pb_all).to_csv(f"{W}/prep/GSE276503_pseudobulk.csv")
json.dump({"dataset":"GSE276503","cells":int(total_cells),
 "design":"29 usable samples (donor structure unstated; each sample = donor unit): control n=2 (HC_IT01/02) vs CRSwNP n=27 (IT/MT/NP) -> composition-only (control arm <3)",
 "upgrade":"pending->core: GSM titles confirm nasal CRS tissue (PMID 39358540)",
 "corrupt_at_source":"CRSIT11, CRSIT14, CRSIT15, CRSMT22: deposited matrix.mtx.gz truncated mid-file (tar member size == filelist size; error 'Invalid integer value'); excluded",
 "file_anomaly":"GSM8499604 files lack donor infix (matrix=GSM8499604_matrix.mtx.gz; barcodes=GSM8499604_CRSNP33_*); donor NP33 from GSM title CRSwNP_NP33",
 "method":"parts clustered separately (memory), labels harmonized by unified marker vocabulary; targeted gene checks resolved ties",
 "note":"basal clusters in parts 0/2 co-express S100A8 (neutrophil-like signal); p1_cl19 CD3D+/CD8A+/LYZhi kept ambiguous (possible doublet)"},
 open(f"{W}/prep/GSE276503_runlog.json","w"), indent=1)
print("OK total cells:", total_cells)
