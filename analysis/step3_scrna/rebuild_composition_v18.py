#!/usr/bin/env python3
"""Round-4 (v1.8) composition rebuild.

The deposited long-format composition_all.csv was assembled with the
'secretory_club' label collapsed into 'secretory' for GSE235711 and
GSE156285 (each donor carried two rows named 'secretory'). The annotated
h5ad objects retain the two labels distinctly. This script rebuilds the
long composition table directly from the objects' obs (donor, group,
celltype), asserts structural integrity, and diffs against the old table:
after summing secretory + secretory_club in the new table, every row must
reproduce the old table exactly (tolerance 1e-9); all other rows must be
identical. Output: composition_all_v1.8.csv (the old file is archived by
the caller, not here).
"""
import h5py, glob, os, sys
import numpy as np, pandas as pd

W = os.path.dirname(os.path.abspath(__file__))

def obs_col(obs, name):
    node = obs[name]
    if isinstance(node, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else c for c in node['categories'][:]]
        codes = node['codes'][:]
        return np.array([cats[c] if c >= 0 else 'NA' for c in codes], dtype=object)
    arr = node[:]
    return np.array([a.decode() if isinstance(a, bytes) else a for a in arr], dtype=object)

RENAME = {'CD4_T_type2': 'CD4_Th2'}  # assembly-time harmonisation rename used by the old table

def composition_from_obs(frames):
    """frames: list of (donor, group, celltype) arrays -> long composition rows
    on the full per-dataset grid (every label of the dataset for every donor,
    fraction 0 where absent), matching the old table's convention."""
    obs = pd.concat([pd.DataFrame({'donor': d, 'group': g, 'celltype': c})
                     for d, g, c in frames], ignore_index=True)
    obs['celltype'] = obs['celltype'].replace(RENAME)
    cnt = obs.groupby(['donor', 'group', 'celltype']).size().rename('n').reset_index()
    tot = cnt.groupby(['donor', 'group'])['n'].transform('sum')
    cnt['fraction'] = cnt['n'] / tot
    donors = obs[['donor', 'group']].drop_duplicates()
    grid = donors.merge(pd.DataFrame({'celltype': sorted(obs['celltype'].unique())}),
                        how='cross')
    out = grid.merge(cnt[['donor', 'celltype', 'fraction']],
                     on=['donor', 'celltype'], how='left')
    out['fraction'] = out['fraction'].fillna(0.0)
    return out[['donor', 'group', 'celltype', 'fraction']]

parts = []
for path in sorted(glob.glob(f'{W}/objects/*_annotated.h5ad')):
    base = os.path.basename(path).replace('_annotated.h5ad', '')
    ds = base.split('_part')[0] if '_part' in base else base
    with h5py.File(path, 'r') as f:
        obs = f['obs']
        frames = [(obs_col(obs, 'donor'), obs_col(obs, 'group'), obs_col(obs, 'celltype'))]
    comp = composition_from_obs(frames)
    comp.insert(0, 'dataset', ds)
    parts.append(comp)
    print(f'{base}: {len(comp)} rows, {comp.donor.nunique()} donors')

# GSE276503 was deposited as four part-objects but is one dataset: re-grid on
# the union label set across parts (the old table's convention).
i276 = [i for i, p in enumerate(parts) if (p['dataset'] == 'GSE276503').all()]
if i276:
    sub = pd.concat([parts[i] for i in i276], ignore_index=True)
    donors = sub[['donor', 'group']].drop_duplicates()
    grid = donors.merge(pd.DataFrame({'celltype': sorted(sub['celltype'].unique())}),
                        how='cross')
    sub = grid.merge(sub[['donor', 'group', 'celltype', 'fraction']],
                     on=['donor', 'group', 'celltype'], how='left')
    sub['fraction'] = sub['fraction'].fillna(0.0)
    sub.insert(0, 'dataset', 'GSE276503')
    parts = [p for i, p in enumerate(parts) if i not in i276] + [sub[['dataset', 'donor', 'group', 'celltype', 'fraction']]]

new = pd.concat(parts, ignore_index=True)[['dataset', 'donor', 'group', 'celltype', 'fraction']]

# ---- structural assertions ----
dup = new.duplicated(subset=['dataset', 'donor', 'celltype']).sum()
assert dup == 0, f'{dup} duplicate (dataset,donor,celltype) rows remain'
sums = new.groupby(['dataset', 'donor'])['fraction'].sum()
assert np.allclose(sums, 1.0, atol=1e-9), f'fraction sums off: {sums[(sums - 1).abs() > 1e-9]}'
gm = new.groupby(['dataset', 'donor'])['group'].nunique()
assert (gm == 1).all(), 'a donor maps to >1 group'

# ---- diff against old table ----
old = pd.read_csv(f'{W}/composition_all.csv')
assert list(old.columns) == list(new.columns), (old.columns, new.columns)
# the old table itself carries duplicate 'secretory' keys; aggregate it the
# same way before comparing (a plain merge would cartesian-join the pairs)
old = old.groupby(['dataset', 'donor', 'group', 'celltype'], as_index=False)['fraction'].sum()
n2 = new.copy()
n2['celltype'] = n2['celltype'].replace({'secretory_club': 'secretory'})
n2 = n2.groupby(['dataset', 'donor', 'group', 'celltype'], as_index=False)['fraction'].sum()
m = old.merge(n2, on=['dataset', 'donor', 'group', 'celltype'], suffixes=('_old', '_new'), how='outer', indicator=True)
assert (m['_merge'] == 'both').all(), m[m['_merge'] != 'both'].head(20)
m['adelta'] = (m['fraction_old'] - m['fraction_new']).abs()
if os.environ.get('DEBUG_DIFF'):
    print(m.sort_values('adelta', ascending=False).head(25).to_string())
    sys.exit(0)
delta = m['adelta'].max()
assert delta < 1e-9, f'max fraction delta after merge-back: {delta}'
print(f'old rows {len(old)} -> new rows {len(new)} (delta {len(new) - len(old)})')
print('merge-back max delta:', delta)
print('new distinct celltype values:', new.celltype.nunique(), sorted(new.celltype.unique()))
for ds, sub in new.groupby('dataset'):
    sc = sub[sub.celltype == 'secretory_club']
    if len(sc):
        print(f'  {ds}: secretory_club rows {len(sc)} across {sc.donor.nunique()} donors')

new.to_csv(f'{W}/composition_all_v1.8.csv', index=False)
print('wrote composition_all_v1.8.csv')
