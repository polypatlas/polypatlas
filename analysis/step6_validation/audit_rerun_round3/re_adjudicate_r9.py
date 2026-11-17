#!/usr/bin/env python3
"""Round-3 re-adjudication: move R9 (study-level claim) out of the gene-level chain.
Input : analysis/step6_validation/reproduction_audit_v1.6.csv (frozen)
Output: analysis/step6_validation/reproduction_audit_v1.7.csv + denominator_chain_v4.md
Only R9 changes: verdict reproduced -> not_testable, note documents the study-level
assessment. Everything else is copied byte-identical per row; the chain is recomputed
from the file itself and asserted against expectations."""
import csv, collections, pathlib

SRC = pathlib.Path('analysis/step6_validation/reproduction_audit_v1.6.csv')
DST = pathlib.Path('analysis/step6_validation/reproduction_audit_v1.7.csv')
rows = list(csv.DictReader(open(SRC, encoding='utf-8')))
fields = list(rows[0].keys())
assert len(rows) == 85, len(rows)

R9_NOTE = ("study-level claim (eosinophilic vs non-eosinophilic transcriptome profiles), "
           "assessed separately by its own criterion: 3 FDR-significant DEGs in the "
           "eosinophilic-vs-non-eosinophilic contrast of GSE72713; excluded from the "
           "gene-level chain because it carries no gene, no expected direction and no "
           "claim-level FDR (round-3 re-adjudication)")

n_changed = 0
for r in rows:
    if r['id'] == 'R9':
        assert r['verdict'] == 'reproduced' and r['gene'] == 'N/A' and r['expected'] == 'N/A'
        assert r['FDR'] == '' and 'observed FDR<0.05 count=3' in r['note']
        old_note = r['note']
        r['verdict'] = 'not_testable'
        r['note'] = (old_note.rstrip('; ') + '; ' + R9_NOTE) if old_note.strip('; ') else R9_NOTE
        n_changed += 1
assert n_changed == 1

with open(DST, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(rows)

# ---- recompute the chain from the written file ----
rows = list(csv.DictReader(open(DST, encoding='utf-8')))
TEST = ('reproduced', 'direction_only', 'not_reproduced', 'not_measured')
DEF = ('reproduced', 'direction_only', 'not_reproduced')
def chain(rs):
    t = [r for r in rs if r['verdict'] in TEST]
    c = collections.Counter(r['verdict'] for r in t)
    dc = c['reproduced'] + c['direction_only']
    meas = [r for r in t if r['verdict'] != 'not_measured']
    return t, c, dc, meas
t, c, dc, meas = chain(rows)
assert len(t) == 60, (len(t), c)
assert (c['reproduced'], c['direction_only'], c['not_reproduced'], c['not_measured']) == (23, 30, 3, 4), c
assert dc == 53
assert len(meas) == 56
mc = collections.Counter(r['verdict'] for r in meas)
assert mc['reproduced'] == 23 and mc['reproduced'] + mc['direction_only'] == 53
ex = [r for r in t if r['match'] == 'exact']
ec = collections.Counter(r['verdict'] for r in ex)
assert len(ex) == 17 and (ec['reproduced'], ec['direction_only'], ec['not_reproduced'], ec['not_measured']) == (9, 7, 0, 1), (len(ex), ec)
px = [r for r in t if r['match'] == 'proxy']
pc = collections.Counter(r['verdict'] for r in px)
assert len(px) == 43 and (pc['reproduced'], pc['direction_only'], pc['not_reproduced'], pc['not_measured']) == (14, 23, 3, 3), (len(px), pc)
gene_def = [r for r in rows if r['verdict'] in DEF and r['gene'] not in ('N/A', '') and 'fraction' not in r['gene']]
gc = collections.Counter(r['verdict'] for r in gene_def)
assert len(gene_def) == 53 and (gc['reproduced'], gc['direction_only'], gc['not_reproduced']) == (22, 28, 3), (len(gene_def), gc)
genes = sorted(set(r['gene'] for r in gene_def)); assert len(genes) == 47, len(genes)
allc = collections.Counter(r['verdict'] for r in rows)
assert allc['not_testable'] == 21 and allc['skipped'] == 4, allc

md = f"""# Denominator chain v4 (release v1.7, round-3 re-adjudication)

Source: `reproduction_audit_v1.7.csv` (85 rows), recomputed by
`audit_rerun_round3/re_adjudicate_r9.py` from the frozen v1.6 file; the only
row changed is R9.

## What changed and why

R9 (GSE72713) is a study-level claim — "eosinophilic and non-eosinophilic
CRSwNP display distinct transcriptome profiles" — with no gene, no expected
direction and no claim-level FDR (gene = N/A, expected = N/A, FDR empty). It
had been graded `reproduced` on the strength of 3 FDR-significant DEGs in the
same contrast. That grade does not meet the §2.6 gene-level definition of
reproduced (correct direction and FDR < 0.05 for the claimed gene), so R9 is
moved out of the gene-level chain: verdict `not_testable`, with the
study-level assessment recorded in its note (3 FDR-significant DEGs in the
eosinophilic-vs-non-eosinophilic contrast). All other rows are unchanged.

## Chain (eligible denominator)

- Total rows: 85
- Testable (eligible): **60** = 23 reproduced / 30 direction_only /
  3 not_reproduced / 4 not_measured
- Not testable: 21 (incl. R9); skipped (pathway-level): 4
- Direction-consistent: **53/60 = 88.3%**

## Measured denominator (excluding the 4 not_measured)

- Measured: **56** = 23 reproduced / 30 direction_only / 3 not_reproduced
- Reproduced: **23/56 = 41.1%**; direction-consistent: **53/56 = 94.6%**

## Strata

- Exact-match: 17 testable = 9 reproduced / 7 direction_only /
  0 not_reproduced / 1 not_measured; the 16 measured exact claims are all
  direction-consistent, none failed to reproduce.
- Proxy: 43 testable = 14 / 23 / 3 / 3; direction-consistent 37/43 = 86.0%.

## Gene-level definitive claims (literature badges)

- 53 claims (22 reproduced / 28 direction_only / 3 not_reproduced) over
  47 genes — unchanged from v1.6 (R9 is not gene-anchored).
- Non-gene definitive claims: 3 (R23-mast reproduced; R32b, R32c
  direction_only).
"""
pathlib.Path('analysis/step6_validation/denominator_chain_v4.md').write_text(md, encoding='utf-8')
print('OK: v1.7 written;', len(rows), 'rows; testable', len(t), dict(c),
      '; exact', len(ex), dict(ec), '; gene-def', len(gene_def), dict(gc))
