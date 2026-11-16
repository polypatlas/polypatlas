#!/usr/bin/env python3
"""Consensus-rule verification + probe-collapse sensitivity consensus.

Step 1: verify the inferred main consensus rule reproduces the frozen
step6 consensus_genes.csv (disease_cohorts_10: 825 up / 196 down) EXACTLY
from the frozen step2 de/ tables, restricted to the evidence-matrix
(protein-coding) gene universe.
Step 2: substitute the 4 microarray disease cohorts' primary contrasts with
the max-variance-collapse results (de_maxvar/) and re-apply the same rule.
"""
import csv, json, statistics, sys
from pathlib import Path

BASE = Path("analysis/step2_bulk")
SENS = BASE / "sensitivity_probe_collapse"
STEP6 = BASE.parent / "step6_validation"

PRIMARY = {  # the 10 disease-vs-control cohorts (same set as f1_meta.R)
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
CHIP = {"GSE36830", "GSE194282", "GSE69093", "GSE70898_70900"}

def load_de(path):
    # a gene counts as MEASURED when logFC is present, even if FDR is NA
    # (DESeq2 padj=NA rows exist in the frozen de/ tables; the original
    #  consensus counted them as measured, verified on boundary genes)
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

# gene universe = evidence matrix rows (protein-coding)
with open(BASE / "evidence_matrix_logFC.csv") as f:
    universe = [r[""] for r in csv.DictReader(f)]
print(f"universe genes: {len(universe)}")

def build(de_dir):
    return {g: load_de(de_dir / f"{PRIMARY[g]}.csv") for g in PRIMARY}

def apply_rule(mats):
    out = {}  # gene -> dict
    for gene in universe:
        vals = [(mats[g][gene][0], mats[g][gene][1]) for g in PRIMARY if gene in mats[g]]
        n_meas = len(vals)
        if n_meas < 6:
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

main_mats = build(BASE / "de")
main_cons = apply_rule(main_mats)
up = {g for g, v in main_cons.items() if v["direction"] == "up"}
down = {g for g, v in main_cons.items() if v["direction"] == "down"}
print(f"recomputed main consensus: up={len(up)} down={len(down)}")

# verify against frozen step6 file
frozen = {}
with open(STEP6 / "consensus_genes.csv") as f:
    for r in csv.DictReader(f):
        if r["version"] == "disease_cohorts_10":
            frozen[r["gene"]] = r
f_up = {g for g, r in frozen.items() if r["direction"] == "up"}
f_down = {g for g, r in frozen.items() if r["direction"] == "down"}
print(f"frozen: up={len(f_up)} down={len(f_down)}")
assert up == f_up and down == f_down, (
    f"RULE MISMATCH: only-recomputed up={sorted(up - f_up)[:10]} down={sorted(down - f_down)[:10]}; "
    f"only-frozen up={sorted(f_up - up)[:10]} down={sorted(f_down - down)[:10]}")
# verify the per-gene count columns too
bad = [g for g in frozen if int(frozen[g]["n_cohorts_measured"]) != main_cons[g]["n_cohorts_measured"]
       or int(frozen[g]["n_consistent"]) != main_cons[g]["n_consistent"]
       or int(frozen[g]["n_sig_FDR005"]) != main_cons[g]["n_sig_FDR005"]
       or abs(float(frozen[g]["median_logFC"]) - main_cons[g]["median_logFC"]) > 0.0011]
assert not bad, f"count-column mismatches: {bad[:10]}"
print("RULE VERIFIED: gene sets AND count columns reproduce frozen consensus_genes.csv exactly")

# sensitivity consensus: chip cohorts from de_maxvar/, RNA-seq unchanged
sens_mats = {g: (load_de(SENS / "de_maxvar" / f"{PRIMARY[g]}.csv") if g in CHIP else main_mats[g])
             for g in PRIMARY}
sens_cons = apply_rule(sens_mats)
s_up = {g for g, v in sens_cons.items() if v["direction"] == "up"}
s_down = {g for g, v in sens_cons.items() if v["direction"] == "down"}
print(f"sensitivity consensus: up={len(s_up)} down={len(s_down)}")

rows = []
allg = sorted(set(main_cons) | set(sens_cons))
for g in allg:
    m, s = main_cons.get(g), sens_cons.get(g)
    rows.append(dict(gene=g,
                     main_direction=m["direction"] if m else "",
                     sens_direction=s["direction"] if s else "",
                     main_n_sig=m["n_sig_FDR005"] if m else "",
                     sens_n_sig=s["n_sig_FDR005"] if s else "",
                     main_frac=m["frac_consistent"] if m else "",
                     sens_frac=s["frac_consistent"] if s else "",
                     status=("both" if m and s else "main_only" if m else "sens_only")))
with open(SENS / "consensus_sensitivity.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

def stats(a, b):
    inter = len(a & b)
    return dict(main=len(a), sens=len(b), retained=inter,
                retention_of_main=round(inter / len(a), 4) if a else None,
                jaccard=round(inter / len(a | b), 4) if (a | b) else None)
summary = {"up": stats(up, s_up), "down": stats(down, s_down),
           "all": stats(up | down, s_up | s_down),
           "flipped_direction": sorted((up & s_down) | (down & s_up))}
print(json.dumps(summary, indent=1))
with open(SENS / "consensus_overlap.json", "w") as f:
    json.dump(summary, f, indent=1)

# anchor genes across all chip primary contrasts (main vs maxvar)
ANCHOR_PRIMARY = dict(PRIMARY)
ANCHOR_PRIMARY.update({"GSE267331": "GSE267331_dupilumab_post_vs_pre",
                       "GSE150846": "GSE150846_NP_ILC2_vs_blood_ILC2"})
arows = []
for g in ["GSE36830", "GSE194282", "GSE69093", "GSE70898_70900", "GSE267331", "GSE150846"]:
    md = load_de(BASE / "de" / f"{ANCHOR_PRIMARY[g]}.csv")
    sd = load_de(SENS / "de_maxvar" / f"{ANCHOR_PRIMARY[g]}.csv")
    for gene in ["CLC", "IL33"]:
        m, s = md.get(gene), sd.get(gene)
        arows.append(dict(cohort=g, contrast=ANCHOR_PRIMARY[g], gene=gene,
                          main_logFC=round(m[0], 4) if m else "", main_FDR=m[1] if m else "",
                          maxvar_logFC=round(s[0], 4) if s else "", maxvar_FDR=s[1] if s else "",
                          main_sig=bool(m and m[1] is not None and m[1] < 0.05), maxvar_sig=bool(s and s[1] is not None and s[1] < 0.05),
                          same_sign=bool(m and s and (m[0] > 0) == (s[0] > 0))))
    for gene in ["CLC", "IL33"]:
        pass
for gene in ["CLC", "IL33"]:
    arows.append(dict(cohort="CONSENSUS", contrast="disease_cohorts_10", gene=gene,
                      main_logFC="", main_FDR="", maxvar_logFC="", maxvar_FDR="",
                      main_sig=("in_main_consensus:" + (main_cons[gene]["direction"] if gene in main_cons else "no")),
                      maxvar_sig=("in_sens_consensus:" + (sens_cons[gene]["direction"] if gene in sens_cons else "no")),
                      same_sign=(gene in main_cons) == (gene in sens_cons)))
with open(SENS / "anchor_genes.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(arows[0].keys())); w.writeheader(); w.writerows(arows)
print("anchors written")
