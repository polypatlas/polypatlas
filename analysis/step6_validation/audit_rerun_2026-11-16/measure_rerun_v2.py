"""Corrected measurements (v2): uniform primary contrast NP (n=15) vs HC-IT (n=2);
each test restricted to donors actually measured in the product (per-part label
sets differ: 'secretory'/'monocyte' labels absent in parts 1-2 annotations,
MMP3 absent from parts 0-1 gene sets). Supersedes measure_rerun.py values where
they differ; see extraction_notes.md."""
import numpy as np, pandas as pd, anndata as ad
from scipy.stats import mannwhitneyu
S3 = 'analysis/step3_scrna'
NP = ['CRSNP31','CRSNP32','CRSNP34','CRSNP35','CRSNP36','CRSNP37','CRSNP38','CRSNP39',
      'CRSNP40','CRSNP41','CRSNP42','CRSNP43','CRSNP44','CRSNP45','NP33']
HC = ['HCIT01','HCIT02']
def test(a, b, pc):
    a = np.array([x for x in a if not np.isnan(x)], float); b = np.array([x for x in b if not np.isnan(x)], float)
    logfc = np.log2((a.mean()+pc)/(b.mean()+pc))
    p = mannwhitneyu(a, b, alternative='two-sided').pvalue
    v = 'reproduced' if (logfc>0 and p<0.05) else ('direction_only' if logfc>0 else 'not_reproduced')
    return logfc, p, v, len(a), len(b)
pb = pd.read_csv(f'{S3}/prep/GSE276503_pseudobulk.csv')
pb['lib'] = pb.drop(columns=['donor','celltype']).sum(axis=1)
def cpm(donor, ct, gene):
    r = pb[(pb.donor==donor) & (pb.celltype==ct)]
    if len(r)==0: return np.nan
    val = r.iloc[0][gene]
    return val / r.iloc[0]['lib'] * 1e6
for ct, genes in [('secretory', ['S100A8','S100A9','IL1RN','AREG','EREG']), ('fibroblast', ['MMP3'])]:
    for g in genes:
        logfc, p, v, na, nb = test([cpm(d, ct, g) for d in NP], [cpm(d, ct, g) for d in HC], pc=0.1)
        print(f'PB {ct} {g}: measured NP n={na}, HC n={nb}  logFC={logfc:.3f} p={p:.4f} verdict={v}')
# monocyte fraction, primary NP vs HC (composition_all.csv)
comp = pd.read_csv(f'{S3}/composition_all.csv')
c = comp[(comp.dataset=='GSE276503') & (comp.celltype=='monocyte')].drop_duplicates(['donor'])
mono = dict(zip(c.donor, c.fraction))
logfc, p, v, na, nb = test([mono[d] for d in NP], [mono[d] for d in HC], pc=1e-4)
print(f'monocyte_fraction NP vs HC: n={na} vs {nb}  logFC={logfc:.3f} p={p:.4f} verdict={v}')
# LY6D+ fraction within secretory, part3 only (only part with secretory label AND LY6D gene AND both arms)
a3 = ad.read_h5ad(f'{S3}/objects/GSE276503_part3_annotated.h5ad')
sec = a3[a3.obs['celltype']=='secretory']
vals = np.asarray(sec[:, 'LY6D'].X.todense()).ravel()
df = pd.DataFrame({'donor': sec.obs['donor'].astype(str).values, 'pos': vals>0})
fr = df.groupby('donor')['pos'].mean()
print('LY6D+ frac within secretory (part3):', {k: round(x,4) for k,x in fr.items()})
logfc, p, v, na, nb = test([fr[d] for d in NP if d in fr], [fr[d] for d in HC if d in fr], pc=1e-4)
print(f'LY6D+ fraction NP vs HC (part3 donors): n={na} vs {nb}  logFC={logfc:.3f} p={p:.4f} verdict={v}')
