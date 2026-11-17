#!/usr/bin/env python3
"""GSE198950-exclusion sensitivity consensus (PolypAtlas supplement).

Verbatim reuse of the rule/loader in
../sensitivity_strict_tissue/strict_consensus.py (itself reusing
../sensitivity_probe_collapse/consensus_sensitivity.py, which verified the
disease-only rule reproduces the frozen step6 consensus 825 up / 196 down
exactly; this script re-asserts that reproduction, plus the strict-tissue
545 up / 225 down at >=5/7, before computing anything new).

Rule (disease-only, unchanged): measured = logFC exists; direction = sign of
median logFC; consistency (same-sign fraction) >= 0.8; same-direction
FDR<0.05 cohorts >= 3. Only the cohort subset and the measured threshold
change.

Cohort sets (primary disease contrasts, see
../sensitivity_strict_tissue/primary_contrasts_inventory.csv):
  MAIN9   = all 10 primary cohorts minus GSE198950 (9 cohorts).
            Frozen primary threshold 6/10 scaled proportionally:
            ceil(0.6 * 9) = 6  ->  measured >= 6/9 (primary);
            variants >= 5/9 and >= 7/9.
  STRICT6 = strict-tissue 7 cohorts minus GSE198950 (6 cohorts:
            GSE136825, GSE338426, GSE179265, GSE36830, GSE194282, GSE72713).
            Strict primary threshold 5/7 scaled proportionally:
            ceil((5/7) * 6) = ceil(4.286) = 5  ->  measured >= 5/6 (primary);
            variants >= 4/6 and >= 6/6.
"""
import csv, json, math, statistics
from pathlib import Path

BASE = Path("analysis/step2_bulk")
OUT = BASE / "sensitivity_gse198950_out"
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
STRICT7 = {g: c for g, c in PRIMARY_ALL.items()
           if g not in {"GSE69093", "GSE313186", "GSE70898_70900"}}
MAIN9 = {g: c for g, c in PRIMARY_ALL.items() if g != "GSE198950"}
STRICT6 = {g: c for g, c in STRICT7.items() if g != "GSE198950"}
assert len(PRIMARY_ALL) == 10 and len(STRICT7) == 7
assert len(MAIN9) == 9 and len(STRICT6) == 6
assert set(STRICT6) == {"GSE136825", "GSE338426", "GSE179265",
                        "GSE36830", "GSE194282", "GSE72713"}


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
assert len(universe) == 19266

ALL_MATS = {g: load_de(BASE / "de" / f"{c}.csv") for g, c in PRIMARY_ALL.items()}


def apply_rule(cohorts, min_meas):
    mats = {g: ALL_MATS[g] for g in cohorts}
    out = {}
    for gene in universe:
        vals = [(mats[g][gene][0], mats[g][gene][1]) for g in cohorts if gene in mats[g]]
        n_meas = len(vals)
        if n_meas < min_meas:
            continue
        med = statistics.median(v[0] for v in vals)
        if med == 0:
            continue
        direction = 1 if med > 0 else -1
        n_cons = sum(1 for lfc, _ in vals if lfc * direction > 0)
        n_sig = sum(1 for lfc, fdr in vals
                    if lfc * direction > 0 and fdr is not None and fdr < 0.05)
        if n_cons / n_meas >= 0.8 and n_sig >= 3:
            out[gene] = dict(direction="up" if direction > 0 else "down",
                             n_cohorts_measured=n_meas, n_consistent=n_cons,
                             frac_consistent=round(n_cons / n_meas, 4),
                             n_sig_FDR005=n_sig, median_logFC=round(med, 3))
    return out


def split(cons):
    return ({g for g, v in cons.items() if v["direction"] == "up"},
            {g for g, v in cons.items() if v["direction"] == "down"})


def anchor_diag(cohorts, gene, cons):
    """CLC/IL33 status in a set: consensus row if member, else raw tallies
    (same diagnostics block as strict_consensus.py)."""
    v = cons.get(gene)
    if v:
        return dict(in_consensus=True, **v)
    vals = [(ALL_MATS[g][gene][0], ALL_MATS[g][gene][1])
            for g in cohorts if gene in ALL_MATS[g]]
    if not vals:
        return dict(in_consensus=False, n_measured=0)
    med = statistics.median(x[0] for x in vals)
    d = 1 if med > 0 else -1
    return dict(in_consensus=False,
                n_measured=len(vals), median_logFC=round(med, 3),
                direction_by_median="up" if d > 0 else "down",
                n_consistent=sum(1 for l, _ in vals if l * d > 0),
                frac_consistent=round(sum(1 for l, _ in vals if l * d > 0) / len(vals), 4),
                n_sig_same_dir=sum(1 for l, fd in vals
                                   if l * d > 0 and fd is not None and fd < 0.05))


def compare(ref_up, ref_down, new_up, new_down, ref_name):
    def stats(a, b):
        inter = len(a & b)
        return dict(reference=len(a), new=len(b), intersection=inter,
                    retention_of_reference=round(inter / len(a), 4) if a else None,
                    jaccard=round(inter / len(a | b), 4) if (a | b) else None)
    flipped = sorted((ref_up & new_down) | (ref_down & new_up))
    return {
        "reference": ref_name,
        "up": stats(ref_up, new_up), "down": stats(ref_down, new_down),
        "all": stats(ref_up | ref_down, new_up | new_down),
        "direction_flipped_genes": flipped,
        "n_direction_flipped": len(flipped),
        "new_only_up": sorted(new_up - ref_up), "new_only_down": sorted(new_down - ref_down),
        "reference_only_up_count": len(ref_up - new_up),
        "reference_only_down_count": len(ref_down - new_down),
    }


# ---- frozen references -----------------------------------------------------
frozen = {}
with open(STEP6 / "consensus_genes.csv") as f:
    for r in csv.DictReader(f):
        if r["version"] == "disease_cohorts_10":
            frozen[r["gene"]] = r
f_up = {g for g, r in frozen.items() if r["direction"] == "up"}
f_down = {g for g, r in frozen.items() if r["direction"] == "down"}
print(f"frozen main consensus: up={len(f_up)} down={len(f_down)}")
assert (len(f_up), len(f_down)) == (825, 196)

strict_frozen = {}
with open(STRICT_DIR / "strict_consensus.csv") as f:
    for r in csv.DictReader(f):
        strict_frozen[r["gene"]] = r
sf_up = {g for g, r in strict_frozen.items() if r["direction"] == "up"}
sf_down = {g for g, r in strict_frozen.items() if r["direction"] == "down"}
print(f"frozen strict consensus: up={len(sf_up)} down={len(sf_down)}")
assert (len(sf_up), len(sf_down)) == (545, 225)

# ---- sanity: verbatim rule reproduces both frozen lists --------------------
chk10 = apply_rule(list(PRIMARY_ALL), 6)
c10_up, c10_down = split(chk10)
assert c10_up == f_up and c10_down == f_down, "10-cohort reproduction FAILED"
chk7 = apply_rule(list(STRICT7), 5)
c7_up, c7_down = split(chk7)
assert c7_up == sf_up and c7_down == sf_down, "strict-7 reproduction FAILED"
print("ASSERTION PASSED: rule reproduces frozen 825/196 (10 cohorts, >=6/10) "
      "and strict 545/225 (7 cohorts, >=5/7) exactly")

# ---- the two new sets --------------------------------------------------------
SETS = {
    "main9": dict(cohorts=MAIN9, primary=6, variants=(5, 6, 7),
                  primary_label=">=6/9 (ceil(6/10*9)=ceil(5.4)=6)"),
    "strict6": dict(cohorts=STRICT6, primary=5, variants=(4, 5, 6),
                    primary_label=">=5/6 (ceil(5/7*6)=ceil(4.286)=5)"),
}
summary = {"threshold_rule": {
    "main9": "measured >=6/9 primary; proportional scaling of frozen 6/10 "
             "(0.6*9=5.4, ceiling 6); variants >=5/9, >=7/9",
    "strict6": "measured >=5/6 primary; proportional scaling of strict >=5/7 "
               "((5/7)*6=4.286, ceiling 5); variants >=4/6, >=6/6",
    "unchanged": "consistency >=0.8; same-direction FDR<0.05 cohorts >=3; "
                 "direction by median logFC; measured = logFC present"},
    "cohorts": {k: sorted(v["cohorts"]) for k, v in SETS.items()},
    "excluded": "GSE198950 (GSE198950_CRSwNP_vs_NonCRS) removed from both sets",
    "sets": {}}

for name, cfg in SETS.items():
    cohorts = list(cfg["cohorts"])
    results = {}
    for thr in cfg["variants"]:
        cons = apply_rule(cohorts, thr)
        up, down = split(cons)
        results[thr] = (cons, up, down)
        print(f"{name} measured >={thr}/{len(cohorts)}: up={len(up)} down={len(down)} "
              f"total={len(cons)}")
    main_cons, m_up, m_down = results[cfg["primary"]]

    # write primary consensus, same column structure as frozen/strict tables
    rows = [dict(direction=v["direction"], gene=g,
                 **{k: v[k] for k in ("n_cohorts_measured", "n_consistent",
                                      "frac_consistent", "n_sig_FDR005", "median_logFC")})
            for g, v in sorted(main_cons.items(),
                               key=lambda kv: (kv[1]["direction"] != "up", kv[0]))]
    fn = "main9_consensus.csv" if name == "main9" else "strict6_consensus.csv"
    with open(OUT / fn, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["direction", "gene", "n_cohorts_measured",
                                          "n_consistent", "frac_consistent",
                                          "n_sig_FDR005", "median_logFC"])
        w.writeheader(); w.writerows(rows)

    entry = {
        "n_cohorts": len(cohorts),
        "primary_threshold": cfg["primary_label"],
        "threshold_variants": {
            f"ge{t}": {"up": len(results[t][1]), "down": len(results[t][2]),
                       "total": len(results[t][0])} for t in cfg["variants"]},
        "anchors_primary": {gene: anchor_diag(cohorts, gene, main_cons)
                            for gene in ("CLC", "IL33")},
        "comparisons": {},
    }
    if name == "main9":
        entry["comparisons"]["vs_frozen_main_10"] = compare(
            f_up, f_down, m_up, m_down, "frozen disease_cohorts_10 (825 up / 196 down)")
    else:
        entry["comparisons"]["vs_frozen_strict_7"] = compare(
            sf_up, sf_down, m_up, m_down, "frozen strict-tissue >=5/7 (545 up / 225 down)")
        entry["comparisons"]["vs_frozen_main_10"] = compare(
            f_up, f_down, m_up, m_down, "frozen disease_cohorts_10 (825 up / 196 down)")
    summary["sets"][name] = entry

def loss_diag(ref_dir_map, cohorts, min_meas):
    """For reference-consensus genes lost in the new set, the first failing
    criterion when re-evaluated in the new cohort subset under the REFERENCE
    direction (measured -> median direction -> consistency -> n_sig)."""
    reasons = {}
    for gene, ref_dir in ref_dir_map.items():
        vals = [(ALL_MATS[g][gene][0], ALL_MATS[g][gene][1])
                for g in cohorts if gene in ALL_MATS[g]]
        d = 1 if ref_dir == "up" else -1
        if len(vals) < min_meas:
            reason = "measured_below_threshold"
        else:
            med = statistics.median(x[0] for x in vals)
            if med == 0 or (med > 0) != (d > 0):
                reason = "median_direction_reversed"
            elif sum(1 for l, _ in vals if l * d > 0) / len(vals) < 0.8:
                reason = "consistency_below_0.8"
            elif sum(1 for l, fd in vals
                     if l * d > 0 and fd is not None and fd < 0.05) < 3:
                reason = "n_sig_below_3"
            else:
                reason = "passes_all_criteria_ERROR"
        reasons.setdefault(reason, []).append(gene)
    return {k: {"count": len(v), "genes": sorted(v)}
            for k, v in sorted(reasons.items())}


# loss diagnostics for the two primary comparisons (lost genes only,
# membership recomputed directly by set arithmetic)
m9 = summary["sets"]["main9"]
main9_cons = apply_rule(list(MAIN9), 6)
m9_up, m9_down = split(main9_cons)
lost_main = {g: "up" for g in f_up - m9_up} | {g: "down" for g in f_down - m9_down}
m9["comparisons"]["vs_frozen_main_10"]["loss_diagnostics"] = \
    loss_diag(lost_main, list(MAIN9), 6)

s6 = summary["sets"]["strict6"]
strict6_cons = apply_rule(list(STRICT6), 5)
s6_up, s6_down = split(strict6_cons)
assert (s6_up, s6_down) == (set(s6["comparisons"]["vs_frozen_strict_7"]["new_only_up"]) | (sf_up & s6_up),
                            set(s6["comparisons"]["vs_frozen_strict_7"]["new_only_down"]) | (sf_down & s6_down))
lost_strict = {g: "up" for g in sf_up - s6_up} | {g: "down" for g in sf_down - s6_down}
s6["comparisons"]["vs_frozen_strict_7"]["loss_diagnostics"] = \
    loss_diag(lost_strict, list(STRICT6), 5)

# anchor per-cohort values (all 10 primary cohorts) for CLC / IL33
with open(OUT / "anchor_cohort_values.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gene", "cohort", "contrast", "logFC", "FDR"])
    for gene in ("CLC", "IL33"):
        for g, c in PRIMARY_ALL.items():
            if gene in ALL_MATS[g]:
                lfc, fdr = ALL_MATS[g][gene]
                w.writerow([gene, g, c, round(lfc, 4),
                            "" if fdr is None else f"{fdr:.6g}"])
            else:
                w.writerow([gene, g, c, "", ""])

print(json.dumps(summary, indent=1, ensure_ascii=False))
with open(OUT / "consensus_comparison.json", "w") as f:
    json.dump(summary, f, indent=1, ensure_ascii=False)
print("wrote main9_consensus.csv + strict6_consensus.csv + consensus_comparison.json "
      "+ anchor_cohort_values.csv")
