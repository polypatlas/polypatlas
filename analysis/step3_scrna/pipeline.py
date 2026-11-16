"""Step 3 common single-cell processing pipeline."""
import os, gc, gzip, json, warnings
import numpy as np, pandas as pd, scipy.io as sio, scipy.sparse as sp
import anndata as ad

MARKERS = {
    # lineage
    "EPI": ["EPCAM","KRT8","KRT18"],
    "IMM": ["PTPRC"],
    "STR": ["COL1A1","DCN","LUM"],
    # epithelial subtypes
    "basal": ["KRT5","TP63","KRT15"],
    "ciliated": ["FOXJ1","TPPP3","DNAH5","C11orf88"],
    "secretory_club": ["SCGB1A1","SCGB3A1"],
    "goblet": ["MUC5AC","MUC2","TFF3","SPDEF"],
    "glandular": ["LTF","LYZ","AZGP1","PIP","DMBT1"],
    "ionocyte": ["FOXI1","CFTR","ASCL3"],
    "tuft": ["POU2F3","TRPM5","LRMP","GNAT3"],
    "deuterosomal": ["DEUP1","FOXN4","CDC20B"],
    "neuroendocrine": ["ASCL1","CALCA","GRP"],
    # immune subtypes
    "CD4_T": ["CD4","IL7R","CD3D"],
    "CD8_T": ["CD8A","CD8B","CD3D"],
    "Treg": ["FOXP3","IL2RA","IKZF2"],
    "Th2": ["GATA3","IL13","IL5","IL4"],
    "T_other": ["CD3D","CD3E","TRAC"],
    "NK": ["NCAM1","NKG7","GNLY","KLRD1"],
    "ILC2": ["IL7R","GATA3","PTGDR2","IL1RL1","KLRB1"],
    "B": ["CD79A","CD79B","MS4A1"],
    "plasma": ["MZB1","SDC1","IGHG1","JCHAIN","XBP1"],
    "mast": ["TPSAB1","CPA3","KIT","HDC","MS4A2"],
    "macrophage": ["CD68","CD163","MARCO","MSR1"],
    "monocyte": ["CD14","FCN1","S100A8"],
    "cDC": ["CD1C","CLEC9A","XCR1","FLT3"],
    "pDC": ["LILRA4","IRF7","CLEC4C"],
    "neutrophil": ["CSF3R","FCGR3B","S100A8","S100A9"],
    "eosinophil": ["CLC","EPX","PRG2","SIGLEC8"],
    # stromal subtypes
    "fibroblast": ["COL1A1","DCN","LUM","COL3A1"],
    "myofibroblast": ["ACTA2","TAGLN","MYL9"],
    "endothelial": ["PECAM1","VWF","CLDN5"],
    "pericyte": ["RGS5","PDGFRB","NOTCH3"],
}
CTRL_ASSERT = {
    "TPI1": "housekeeping TPI1 detectable",
    "GAPDH": "housekeeping GAPDH detectable",
    "EPCAM": "epithelial marker EPCAM",
    "KRT8": "epithelial marker KRT8",
    "PTPRC": "immune marker PTPRC",
}

def qc_adata(a, min_genes=200, max_genes=6000, max_mt=15.0, min_cells=3):
    a.var_names_make_unique()
    a.obs_names_make_unique()
    mt = np.asarray(a.var_names.str.upper().str.startswith("MT-"))
    X = a.X.tocsr() if sp.issparse(a.X) else a.X
    tot = np.asarray(X.sum(1)).ravel()
    a.obs["pct_mt"] = (np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(tot, 1) * 100) if mt.any() else 0.0
    a.obs["n_genes"] = np.asarray((X > 0).sum(1)).ravel()
    a.obs["n_counts"] = tot
    keep = (a.obs.n_genes >= min_genes) & (a.obs.n_genes <= max_genes) & (a.obs.pct_mt <= max_mt)
    a = a[keep].copy()
    keep_g = (a.X > 0).sum(0).A.ravel() >= min_cells if sp.issparse(a.X) else (a.X > 0).sum(0) >= min_cells
    a = a[:, keep_g].copy()
    return a

def positive_controls(a):
    """Return dict gene->fraction_of_cells_detected; print FAIL if any core gene is all-zero."""
    res, fails = {}, []
    X = a.X.tocsr() if sp.issparse(a.X) else a.X
    for g, desc in CTRL_ASSERT.items():
        if g in a.var_names:
            frac = float((X[:, a.var_names.get_loc(g)] > 0).mean())
            res[g] = round(frac, 4)
            if frac == 0:
                fails.append(g)
        else:
            res[g] = "ABSENT_FROM_VAR"
            fails.append(g + "(absent)")
    return res, fails

def norm_log(a):
    import scanpy as sc
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a

def cluster(a, n_pcs=40, res=0.8, seed=0):
    import scanpy as sc
    sc.pp.highly_variable_genes(a, n_top_genes=2000, subset=False, flavor="seurat")
    # keep lineage markers even if not highly variable (sorted preps: uniform markers drop out of HVG)
    keep_markers = a.var_names.isin({m for genes in MARKERS.values() for m in genes})
    sel = a.var["highly_variable"].values | np.asarray(keep_markers)
    a = a[:, sel].copy()
    sc.pp.pca(a, n_comps=min(n_pcs, a.n_vars - 1, a.n_obs - 1), svd_solver="arpack")
    sc.pp.neighbors(a, n_pcs=min(n_pcs, a.obsm["X_pca"].shape[1] - 1))
    sc.tl.leiden(a, resolution=res, random_state=seed)
    return a

def signature_scores(a):
    """Mean log1p expression of each marker signature per leiden cluster -> DataFrame.
    Sparse-safe: never densifies the full matrix."""
    import pandas as pd
    cl = a.obs["leiden"].astype(str).values
    X = a.X.tocsr() if sp.issparse(a.X) else a.X
    uniq = pd.unique(cl)
    out = {}
    for sig, genes in MARKERS.items():
        g = [x for x in genes if x in a.var_names]
        if not g:
            continue
        idx = a.var_names.get_indexer(g)
        sub = X[:, idx]
        vals = {}
        for c in uniq:
            mask = cl == c
            n = int(mask.sum())
            vals[c] = float(np.asarray(sub[mask].sum()).ravel().sum() / (n * len(g))) if n else 0.0
        out[sig] = pd.Series(vals)
    return pd.DataFrame(out).fillna(0)

def save_h5ad(a, path):
    a.write_h5ad(path, compression="gzip")

def composition_table(a, donor_col="donor", group_col="group", label_col="celltype"):
    ct = pd.crosstab(a.obs[donor_col], a.obs[label_col], normalize="index")
    meta = a.obs[[donor_col, group_col]].drop_duplicates().set_index(donor_col)
    ct = ct.join(meta)
    ct.index.name = "donor"
    return ct.reset_index()

def pseudobulk(a, donor_col="donor", label_col="celltype"):
    """Donor x celltype pseudobulk counts from raw counts stored in a.raw or a.layers['counts']."""
    src = a.raw.to_adata() if a.raw is not None else a
    X = src.X.tocsr() if sp.issparse(src.X) else src.X
    genes = list(src.var_names)
    keys = list(zip(a.obs[donor_col].astype(str), a.obs[label_col].astype(str)))
    df_idx = pd.MultiIndex.from_tuples(sorted(set(keys)), names=["donor", "celltype"])
    mat = np.zeros((len(df_idx), len(genes)), dtype=np.float64)
    pos = {k: i for i, k in enumerate(df_idx)}
    obs_keys = a.obs[[donor_col, label_col]].astype(str)
    for j in range(X.shape[0]):
        k = (obs_keys.iloc[j, 0], obs_keys.iloc[j, 1])
        mat[pos[k]] += np.asarray(X[j].todense()).ravel() if sp.issparse(X) else X[j]
    return pd.DataFrame(mat, index=df_idx, columns=genes)
