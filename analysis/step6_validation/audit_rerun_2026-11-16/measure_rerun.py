"""Rerun measurements for audit fix round C (2026-11-16).
GSE276503 (PMID 39438439) claims measured against step3 products:
 - composition: analysis/step3_scrna/composition_all.csv
 - pseudobulk counts: analysis/step3_scrna/prep/GSE276503_pseudobulk.csv
 - subset fractions: analysis/step3_scrna/objects/GSE276503_part{0-3}_annotated.h5ad
Verdict rule mirrors frozen protocol: direction correct & p<0.05 -> reproduced;
direction correct & p>=0.05 -> direction_only; wrong direction -> not_reproduced;
gene absent from product -> not_measured. p = two-sided Mann-Whitney U (donor level),
stored in the FDR column per R23-mast precedent (raw donor-level p).
Primary contrast pre-chosen: NP (n=15) vs HC-IT (n=2), the paper's disease-vs-healthy read.
"""
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import mannwhitneyu

S3 = 'analysis/step3_scrna'
NP = ['CRSNP31','CRSNP32','CRSNP34','CRSNP35','CRSNP36','CRSNP37','CRSNP38','CRSNP39',
      'CRSNP40','CRSNP41','CRSNP42','CRSNP43','CRSNP44','CRSNP45','NP33']
HC = ['HCIT01','HCIT02']
CRSIT = ['CRSIT12','CRSIT13','CRSIT16','CRSIT17','CRSIT18','CRSIT19']

def test(vals_np, vals_hc, pc):
    a, b = np.asarray(vals_np, float), np.asarray(vals_hc, float)
    logfc = np.log2((a.mean()+pc)/(b.mean()+pc))
    p = mannwhitneyu(a, b, alternative='two-sided').pvalue if a.size>0 and b.size>0 else np.nan
    verdict = 'reproduced' if (logfc>0 and p<0.05) else ('direction_only' if logfc>0 else 'not_reproduced')
    return logfc, p, verdict

out = []
# ---- composition: monocyte fraction (claim: monocytes greater in NP than other tissues)
comp = pd.read_csv(f'{S3}/composition_all.csv')
c = comp[(comp.dataset=='GSE276503') & (comp.celltype=='monocyte')].drop_duplicates(['donor'])
mono = dict(zip(c.donor, c.fraction))
others = [d for d in mono if d not in NP]
logfc, p, v = test([mono[d] for d in NP], [mono[d] for d in others], pc=1e-4)
out.append(('monocyte_fraction | NP vs all-other-tissues', [round(mono[d],4) for d in NP], [round(mono[d],4) for d in others], logfc, p, v))
logfc2, p2, v2 = test([mono[d] for d in NP], [mono[d] for d in HC], pc=1e-4)
out.append(('monocyte_fraction | NP vs HC only (secondary)', [round(mono[d],4) for d in NP], [round(mono[d],4) for d in HC], logfc2, p2, v2))

# ---- pseudobulk expression (CPM within donor x celltype)
pb = pd.read_csv(f'{S3}/prep/GSE276503_pseudobulk.csv')
pb['lib'] = pb.drop(columns=['donor','celltype']).sum(axis=1)
def cpm(donor, celltype, gene):
    r = pb[(pb.donor==donor) & (pb.celltype==celltype)]
    if len(r)==0 or gene not in pb.columns: return np.nan
    return r.iloc[0][gene] / r.iloc[0]['lib'] * 1e6
def cpm_multi(donor, celltypes, gene):
    r = pb[(pb.donor==donor) & (pb.celltype.isin(celltypes))]
    if len(r)==0: return np.nan
    return r[gene].sum() / r['lib'].sum() * 1e6

plan = [('secretory', ['S100A8','S100A9','IL1RN','AREG','EREG']),
        ('fibroblast', ['CXCL1','CXCL2','CXCL3','CXCL5','CXCL6','CXCL8','CCL5','CCL7','CCL8','CCL11','MMP3','LIF'])]
for ct, genes in plan:
    for g in genes:
        if g not in pb.columns:
            out.append((f'PB {ct} {g}', 'ABSENT', 'ABSENT', np.nan, np.nan, 'not_measured')); continue
        a = [cpm(d, ct, g) for d in NP]; b = [cpm(d, ct, g) for d in HC]
        logfc, p, v = test(a, b, pc=0.1)
        out.append((f'PB {ct} {g}', [round(x,2) for x in a], [round(x,2) for x in b], logfc, p, v))
# IL1B in combined myeloid (monocyte+macrophage+pDC)
a = [cpm_multi(d, ['monocyte','macrophage','pDC'], 'IL1B') for d in NP]
b = [cpm_multi(d, ['monocyte','macrophage','pDC'], 'IL1B') for d in HC]
logfc, p, v = test(a, b, pc=0.1)
out.append(('PB myeloid(mon+mac+pDC) IL1B', [round(x,2) for x in a], [round(x,2) for x in b], logfc, p, v))
# descriptive: S100A8/A9 mean CPM by epithelial class (main-source check)
epi = ['secretory','basal','ciliated','goblet','glandular','ionocyte']
for g in ['S100A8','S100A9']:
    means = {ct: np.nanmean([cpm(d, ct, g) for d in NP+HC+CRSIT]) for ct in epi}
    out.append((f'DESCRIPTIVE mean CPM {g} by epithelial class', means, '', np.nan, np.nan, ''))

# ---- subset fractions from annotated objects (detection: count>0 within parent class)
def subset_fraction(gene, parent, parts):
    fr = {}
    for p in parts:
        a = ad.read_h5ad(f'{S3}/objects/GSE276503_part{p}_annotated.h5ad')
        if gene not in a.var_names: continue
        x = a[a.obs['celltype']==parent]
        if x.shape[0]==0: continue
        vals = np.asarray(x[:, gene].X.todense()).ravel() if hasattr(x[:, gene].X, 'todense') else np.asarray(x[:, gene].X).ravel()
        df = pd.DataFrame({'donor': x.obs['donor'].astype(str).values, 'pos': vals>0})
        for d, grp in df.groupby('donor'):
            n, k = len(grp), grp['pos'].sum()
            if d in fr: fr[d] = (fr[d][0]+k, fr[d][1]+n)
            else: fr[d] = (k, n)
        del a, x
    return {d: k/n for d, (k, n) in fr.items()}
ly = subset_fraction('LY6D', 'secretory', [1,2,3])   # part0 lacks LY6D; NP+HC donors all in parts1-3
a = [ly[d] for d in NP if d in ly]; b = [ly[d] for d in HC if d in ly]
logfc, p, v = test(a, b, pc=1e-4)
out.append(('LY6D+ fraction within secretory | NP vs HC (parts1-3)', [round(x,4) for x in a], [round(x,4) for x in b], logfc, p, v))
ido = subset_fraction('IDO1', 'fibroblast', [0,1,2,3])
a = [ido[d] for d in NP if d in ido]; b = [ido[d] for d in HC if d in ido]
logfc, p, v = test(a, b, pc=1e-4)
out.append(('IDO1+ fraction within fibroblast | NP vs HC (parts0-3)', [round(x,4) for x in a], [round(x,4) for x in b], logfc, p, v))
# secondary: IDO1+ fraction NP vs CRS-IT
a2 = [ido[d] for d in NP if d in ido]; b2 = [ido[d] for d in CRSIT if d in ido]
logfc, p, v = test(a2, b2, pc=1e-4)
out.append(('IDO1+ fraction within fibroblast | NP vs CRS-IT (secondary)', [round(x,4) for x in a2], [round(x,4) for x in b2], logfc, p, v))

for name, va, vb, logfc, p, v in out:
    print(f'{name}\n  NP={va}\n  ref={vb}\n  logFC={logfc if isinstance(logfc,float) else logfc} p={p} verdict={v}')
