#!/usr/bin/env python3
"""GSE36830 donor sensitivity: consensus + DE comparison + marker panel.

Scenarios (disease-only consensus rule, identical to frozen / strict scripts:
measured = logFC present; direction = median logFC sign; consistency >=0.8;
same-direction FDR<0.05 cohorts >=3; measured threshold scaled proportionally
from frozen 6/10, ceiling):
  A_main9    : main 10 cohorts minus GSE36830          -> min measured >=6/9
  A_strict6  : strict 7 cohorts minus GSE36830         -> min measured >=4/6
  B_main10   : main 10, GSE36830 DE replaced by polyp-only DE -> >=6/10 (unchanged)
Comparisons: A_main9 vs frozen main (825/196), A_strict6 vs frozen strict
(545/225), B_main10 vs frozen main. Plus (b) GSE36830 old-vs-new DE key genes
and DEG counts, and (c) marker panel direction recovery without GSE36830.
"""
import csv, json, math, re, statistics
from pathlib import Path

BASE = Path("analysis/step2_bulk")
OUT = BASE / "sensitivity_gse36830_donor"
STEP6 = BASE.parent / "step6_validation"
STRICT_DIR = BASE / "sensitivity_strict_tissue"

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
STRICT7 = [g for g in PRIMARY_ALL if g not in {"GSE69093", "GSE313186", "GSE70898_70900"}]
MARKERS = ["CLC", "CCL26", "CCL13", "IL33", "TSLP", "POSTN", "SERPINB2", "ALOX15"]

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

def load_de_full(path):
    d = {}
    with open(path) as f:
        for r in csv.DictReader(f):
            try:
                d[r["symbol"]] = dict(logFC=float(r["logFC"]), t=float(r["t"]),
                                      p=float(r["p"]), FDR=float(r["FDR"]))
            except (TypeError, ValueError):
                continue
    return d

with open(BASE / "evidence_matrix_logFC.csv") as f:
    universe = [r[""] for r in csv.DictReader(f)]
print(f"universe genes: {len(universe)}")

frozen_mats = {g: load_de(BASE / "de" / f"{c}.csv") for g, c in PRIMARY_ALL.items()}
new36830 = load_de(OUT / "GSE36830_CRSwNP_polypOnly_vs_UTctrl.csv")

def scaled_min(n_cohorts):
    return math.ceil(6 / 10 * n_cohorts)

def apply_rule(mats, cohort_list, min_meas):
    out = {}
    for gene in universe:
        vals = [(mats[g][gene][0], mats[g][gene][1]) for g in cohort_list if gene in mats[g]]
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

def anchor_diag(mats, cohort_list, gene):
    vals = [(mats[g][gene][0], mats[g][gene][1]) for g in cohort_list if gene in mats[g]]
    if not vals:
        return dict(n_measured=0)
    med = statistics.median(x[0] for x in vals)
    d = 1 if med > 0 else -1
    return dict(n_measured=len(vals), median_logFC=round(med, 3),
                n_consistent=sum(1 for l, _ in vals if l * d > 0),
                frac_consistent=round(sum(1 for l, _ in vals if l * d > 0) / len(vals), 4),
                n_sig_same_dir=sum(1 for l, fd in vals if l * d > 0 and fd is not None and fd < 0.05))

# frozen reference sets
frozen_main = {}
with open(STEP6 / "consensus_genes.csv") as f:
    for r in csv.DictReader(f):
        if r["version"] == "disease_cohorts_10":
            frozen_main[r["gene"]] = r["direction"]
frozen_strict = {}
with open(STRICT_DIR / "strict_consensus.csv") as f:
    for r in csv.DictReader(f):
        frozen_strict[r["gene"]] = r["direction"]
assert len([g for g, d in frozen_main.items() if d == "up"]) == 825
assert len([g for g, d in frozen_main.items() if d == "down"]) == 196
assert len([g for g, d in frozen_strict.items() if d == "up"]) == 545
assert len([g for g, d in frozen_strict.items() if d == "down"]) == 225
print("frozen references verified: main 825/196, strict 545/225")

def compare(sens, frozen):
    s_up = {g for g, v in sens.items() if v["direction"] == "up"}
    s_down = {g for g, v in sens.items() if v["direction"] == "down"}
    f_up = {g for g, d in frozen.items() if d == "up"}
    f_down = {g for g, d in frozen.items() if d == "down"}
    def stats(a, b):
        inter = len(a & b)
        return dict(frozen=len(a), sensitivity=len(b), intersection=inter,
                    retention_of_frozen=round(inter / len(a), 4) if a else None,
                    jaccard=round(inter / len(a | b), 4) if (a | b) else None,
                    sensitivity_only=len(b - a), frozen_only=len(a - b))
    flipped = sorted((f_up & s_down) | (f_down & s_up))
    return dict(up=stats(f_up, s_up), down=stats(f_down, s_down),
                all=stats(f_up | f_down, s_up | s_down),
                direction_flipped_genes=flipped), s_up, s_down

scenarios = {}
# A_main9
cohorts = [g for g in PRIMARY_ALL if g != "GSE36830"]
scenarios["A_main9_no36830"] = dict(cohorts=cohorts, mats=frozen_mats, frozen=frozen_main)
# A_strict6
cohorts = [g for g in STRICT7 if g != "GSE36830"]
scenarios["A_strict6_no36830"] = dict(cohorts=cohorts, mats=frozen_mats, frozen=frozen_strict)
# B_main10 polyp-only
mats_b = dict(frozen_mats); mats_b["GSE36830"] = new36830
scenarios["B_main10_polypOnly"] = dict(cohorts=list(PRIMARY_ALL), mats=mats_b, frozen=frozen_main)

summary = {"threshold_rule": "min_measured = ceil(0.6 * n_cohorts); consistency>=0.8; n_sig_same_dir>=3",
           "scenarios": {}}
for name, sc in scenarios.items():
    n = len(sc["cohorts"]); mm = scaled_min(n)
    cons = apply_rule(sc["mats"], sc["cohorts"], mm)
    comp, s_up, s_down = compare(cons, sc["frozen"])
    # threshold variants for transparency
    variants = {}
    for delta in (-1, 0, 1):
        t = mm + delta
        if t < 1 or t > n:
            continue
        c2 = apply_rule(sc["mats"], sc["cohorts"], t)
        variants[f"ge{t}"] = dict(up=sum(1 for v in c2.values() if v["direction"] == "up"),
                                  down=sum(1 for v in c2.values() if v["direction"] == "down"))
    anchors = {}
    for gene in ("CLC", "IL33"):
        if gene in cons:
            anchors[gene] = dict(in_consensus=True, **cons[gene])
        else:
            anchors[gene] = dict(in_consensus=False,
                                 in_frozen_consensus=(gene in sc["frozen"]),
                                 frozen_direction=sc["frozen"].get(gene),
                                 diagnostics=anchor_diag(sc["mats"], sc["cohorts"], gene))
    summary["scenarios"][name] = dict(n_cohorts=n, cohorts=sc["cohorts"], min_measured=mm,
                                      up=len(s_up), down=len(s_down), total=len(cons),
                                      threshold_variants=variants, comparison=comp, anchors=anchors)
    print(f"{name}: n={n} min_meas={mm} up={len(s_up)} down={len(s_down)} flipped={len(comp['direction_flipped_genes'])}")
    # write consensus table
    rows = [dict(direction=v["direction"], gene=g, **{k: v[k] for k in
            ("n_cohorts_measured", "n_consistent", "frac_consistent", "n_sig_FDR005", "median_logFC")})
            for g, v in sorted(cons.items(), key=lambda kv: (kv[1]["direction"] != "up", kv[0]))]
    with open(OUT / f"consensus_{name}.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["direction", "gene", "n_cohorts_measured", "n_consistent",
                                          "frac_consistent", "n_sig_FDR005", "median_logFC"])
        w.writeheader(); w.writerows(rows)

with open(OUT / "consensus_comparison.json", "w") as f:
    json.dump(summary, f, indent=1, ensure_ascii=False)

## ---- (b) GSE36830 old vs new DE: DEG counts + key genes ----
old_full = load_de_full(BASE / "de" / "GSE36830_CRSwNP_vs_UTctrl.csv")
new_full = load_de_full(OUT / "GSE36830_CRSwNP_polypOnly_vs_UTctrl.csv")
deg = dict(old_n_genes=len(old_full), new_n_genes=len(new_full),
           old_deg_FDR005=sum(1 for v in old_full.values() if v["FDR"] is not None and v["FDR"] < 0.05),
           new_deg_FDR005=sum(1 for v in new_full.values() if v["FDR"] is not None and v["FDR"] < 0.05),
           old_deg_up=sum(1 for v in old_full.values() if v["FDR"] is not None and v["FDR"] < 0.05 and v["logFC"] > 0),
           new_deg_up=sum(1 for v in new_full.values() if v["FDR"] is not None and v["FDR"] < 0.05 and v["logFC"] > 0),
           old_deg_down=sum(1 for v in old_full.values() if v["FDR"] is not None and v["FDR"] < 0.05 and v["logFC"] < 0),
           new_deg_down=sum(1 for v in new_full.values() if v["FDR"] is not None and v["FDR"] < 0.05 and v["logFC"] < 0))
# genome-wide logFC correlation old vs new (shared genes)
shared = [g for g in old_full if g in new_full]
import math as _m
xs = [old_full[g]["logFC"] for g in shared]; ys = [new_full[g]["logFC"] for g in shared]
mx, my = statistics.mean(xs), statistics.mean(ys)
r = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / _m.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
deg["shared_genes"] = len(shared); deg["logFC_pearson_r_old_vs_new"] = round(r, 4)
print("DEG comparison:", deg)

key_rows = []
for gene in MARKERS:
    o, nw = old_full.get(gene), new_full.get(gene)
    key_rows.append(dict(gene=gene,
                         old_logFC=round(o["logFC"], 4) if o else "", old_FDR=o["FDR"] if o else "",
                         new_logFC=round(nw["logFC"], 4) if nw else "", new_FDR=nw["FDR"] if nw else "",
                         same_sign=bool(o and nw and (o["logFC"] > 0) == (nw["logFC"] > 0))))
with open(OUT / "gse36830_de_key_genes.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(key_rows[0].keys())); w.writeheader(); w.writerows(key_rows)
with open(OUT / "gse36830_de_deg_counts.json", "w") as f:
    json.dump(deg, f, indent=1)

## ---- (c) marker panel without GSE36830 ----
pat = re.compile(r"([A-Za-z0-9_]+):logFC=(-?[\d.]+),FDR=([\d.eE+-]+)")
panel_rows = []
tot_old = tot_new = 0
with open(STEP6 / "marker_recovery.csv") as f:
    for row in csv.DictReader(f):
        if row["layer"] != "bulk" or row["contrast_class"] != "disease":
            continue
        gene = row["gene"]; expected = 1 if row["expected_direction"] == "up" else -1
        entries = {g: (float(l), float(fd)) for g, l, fd in pat.findall(row["details"])}
        assert len(entries) == 7, (gene, len(entries))
        correct_old = {g: (l * expected > 0) for g, (l, fd) in entries.items()}
        n_old = sum(correct_old.values())
        rest = {g: c for g, c in correct_old.items() if g != "GSE36830"}
        n_new = sum(rest.values())
        tot_old += n_old; tot_new += n_new
        panel_rows.append(dict(gene=gene, expected_direction=row["expected_direction"],
                               frozen_n_correct=n_old, frozen_n_measured=7,
                               gse36830_logFC=entries["GSE36830"][0], gse36830_FDR=entries["GSE36830"][1],
                               gse36830_correct=correct_old["GSE36830"],
                               no36830_n_correct=n_new, no36830_n_measured=6))
assert tot_old == 49, tot_old
panel_rows.append(dict(gene="TOTAL", expected_direction="", frozen_n_correct=tot_old, frozen_n_measured=56,
                       gse36830_logFC="", gse36830_FDR="", gse36830_correct="",
                       no36830_n_correct=tot_new, no36830_n_measured=48))
with open(OUT / "marker_panel_no36830.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(panel_rows[0].keys())); w.writeheader(); w.writerows(panel_rows)
print(f"marker panel: frozen {tot_old}/56 -> without GSE36830 {tot_new}/48")
print("wrote consensus_*.csv, consensus_comparison.json, gse36830_de_key_genes.csv, "
      "gse36830_de_deg_counts.json, marker_panel_no36830.csv")
