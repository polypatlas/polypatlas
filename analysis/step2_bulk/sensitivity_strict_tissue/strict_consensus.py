#!/usr/bin/env python3
"""Strict-tissue sensitivity consensus (supplement B).

Reuses load_de / apply_rule logic from
../sensitivity_probe_collapse/consensus_sensitivity.py (which verified the
disease-only rule reproduces frozen step6 consensus 825 up / 196 down exactly).

Strict set: drop GSE69093 (primary cultured cells), GSE313186 (toxin/deployment
exposure design), GSE70898_70900 (sorted/explant T cells) -> 7 whole-tissue
cohorts. Primary threshold: measured >=5/7 (proportional scaling of 6/10),
consistency >=0.8, same-direction FDR<0.05 cohorts >=3 (unchanged).
Variants reported: measured >=4/7 and >=6/7.
"""
import csv, json, statistics
from pathlib import Path

BASE = Path("analysis/step2_bulk")
OUT = BASE / "sensitivity_strict_tissue"
STEP6 = BASE.parent / "step6_validation"

PRIMARY_ALL = {
    "GSE136825": "GSE136825_polyp_vs_ctrl",
    "GSE338426": "GSE338426_CRSwNP_vs_ctrl",
    "GSE179265": "GSE179265_CRSwNP_vs_ctrl",
    "GSE36830": "GSE36830_CRSwNP_vs_UTctrl",
    "GSE194282": "GSE194282_NP_vs_ctrl",
    "GSE69093": "GSE69093_case_vs_ctrl_celltypeAdj",
    "GSE313186": "GSE313186_exposed_vs_ctrl",
    "GSE70898_70900": "GSE70898_70900_polypT_vs_nasalT",
    "GSE198950": "GSE198950_CRSwNP_vs_NonCRS",
    "GSE72713": "GSE72713_CRSwNP_vs_ctrl",
}
EXCLUDED = {"GSE69093", "GSE313186", "GSE70898_70900"}
STRICT = {g: c for g, c in PRIMARY_ALL.items() if g not in EXCLUDED}
assert len(STRICT) == 7

def load_de(path):
    d = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                lfc = float(r["logFC"])
            except (TypeError, ValueError):
                continue
            try:
                fdr = float(r["FDR"])
            except (TypeError, ValueError):
                fdr = None
            d[r["symbol"]] = (lfc, fdr)
    return d

with open(BASE / "evidence_matrix_logFC.csv") as f:
    universe = [r[""] for r in csv.DictReader(f)]
print(f"universe genes: {len(universe)}")

mats = {g: load_de(BASE / "de" / f"{c}.csv") for g, c in STRICT.items()}

def apply_rule(min_meas):
    out = {}
    for gene in universe:
        vals = [(mats[g][gene][0], mats[g][gene][1]) for g in STRICT if gene in mats[g]]
        n_meas = len(vals)
        if n_meas < min_meas:
            continue
        med = statistics.median(v[0] for v in vals)
        if med == 0:
            continue
        direction = 1 if med > 0 else -1
        n_cons = sum(1 for lfc, _ in vals if lfc * direction > 0)
        n_sig = sum(1 for lfc, fdr in vals if lfc * direction > 0 and fdr is not None and fdr < 0.05)
        if n_cons / n_meas >= 0.8 and n_sig >= 3:
            out[gene] = dict(direction="up" if direction > 0 else "down",
                             n_cohorts_measured=n_meas, n_consistent=n_cons,
                             frac_consistent=round(n_cons / n_meas, 4),
                             n_sig_FDR005=n_sig, median_logFC=round(med, 3))
    return out

results = {}
for thr in (4, 5, 6):
    cons = apply_rule(thr)
    up = {g for g, v in cons.items() if v["direction"] == "up"}
    down = {g for g, v in cons.items() if v["direction"] == "down"}
    results[thr] = (cons, up, down)
    print(f"measured >={thr}/7: up={len(up)} down={len(down)} total={len(cons)}")

main_cons, s_up, s_down = results[5]

# write main strict consensus, same column structure as frozen table
rows = [dict(direction=v["direction"], gene=g, **{k: v[k] for k in
        ("n_cohorts_measured", "n_consistent", "frac_consistent", "n_sig_FDR005", "median_logFC")})
        for g, v in sorted(main_cons.items(), key=lambda kv: (kv[1]["direction"] != "up", kv[0]))]
with open(OUT / "strict_consensus.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["direction", "gene", "n_cohorts_measured", "n_consistent",
                                      "frac_consistent", "n_sig_FDR005", "median_logFC"])
    w.writeheader(); w.writerows(rows)

# frozen comparison
frozen = {}
with open(STEP6 / "consensus_genes.csv") as f:
    for r in csv.DictReader(f):
        if r["version"] == "disease_cohorts_10":
            frozen[r["gene"]] = r
f_up = {g for g, r in frozen.items() if r["direction"] == "up"}
f_down = {g for g, r in frozen.items() if r["direction"] == "down"}
print(f"frozen: up={len(f_up)} down={len(f_down)}")

def stats(a, b):
    inter = len(a & b)
    return dict(frozen=len(a), strict=len(b), intersection=inter,
                retention_of_frozen=round(inter / len(a), 4) if a else None,
                jaccard=round(inter / len(a | b), 4) if (a | b) else None)

flipped = sorted((f_up & s_down) | (f_down & s_up))
summary = {
    "threshold_variants": {f"ge{t}": {"up": len(results[t][1]), "down": len(results[t][2]),
                                      "total": len(results[t][0])} for t in (4, 5, 6)},
    "main_ge5": {"up": stats(f_up, s_up), "down": stats(f_down, s_down),
                 "all": stats(f_up | f_down, s_up | s_down)},
    "direction_flipped_genes": flipped,
    "strict_only_up": sorted(s_up - f_up), "strict_only_down": sorted(s_down - f_down),
    "frozen_only_up_count": len(f_up - s_up), "frozen_only_down_count": len(f_down - s_down),
    "anchors": {},
}
for gene in ("CLC", "IL33"):
    v = main_cons.get(gene)
    summary["anchors"][gene] = (dict(in_strict_consensus=True, **v) if v else
                                dict(in_strict_consensus=False,
                                     in_frozen_consensus=(gene in frozen),
                                     frozen_direction=frozen[gene]["direction"] if gene in frozen else None))
    # near-miss diagnostics for non-members: raw tallies across the 7 cohorts
    if not v:
        vals = [(mats[g][gene][0], mats[g][gene][1]) for g in STRICT if gene in mats[g]]
        if vals:
            med = statistics.median(x[0] for x in vals)
            d = 1 if med > 0 else -1
            summary["anchors"][gene]["diagnostics"] = dict(
                n_measured=len(vals), median_logFC=round(med, 3),
                n_consistent=sum(1 for l, _ in vals if l * d > 0),
                frac_consistent=round(sum(1 for l, _ in vals if l * d > 0) / len(vals), 4),
                n_sig_same_dir=sum(1 for l, fd in vals if l * d > 0 and fd is not None and fd < 0.05))
print(json.dumps(summary, indent=1, ensure_ascii=False))
with open(OUT / "consensus_comparison.json", "w") as f:
    json.dump(summary, f, indent=1, ensure_ascii=False)
print("wrote strict_consensus.csv + consensus_comparison.json")
