#!/usr/bin/env python3
"""F10: treatment-response content pack (existing results only, no new DE)."""
import csv, json, os

BASE = os.path.expanduser("~/workspace/sci-crswnp-resource")
OUT = os.path.join(BASE, "analysis/step9_features")

# (a) marker reversal: 8 marker genes, read straight from the DE table
markers = ["CLC", "CCL26", "CCL13", "IL33", "TSLP", "POSTN", "SERPINB2", "ALOX15"]
de_path = os.path.join(BASE, "analysis/step2_bulk/de/GSE267331_dupilumab_post_vs_pre.csv")
de = {}
with open(de_path) as f:
    for row in csv.DictReader(f):
        g = row["symbol"].strip('"')
        if g in markers:
            de[g] = {k: float(row[k]) if k in ("logFC", "t", "p", "FDR") else row[k]
                     for k in ("logFC", "t", "p", "FDR")}
marker_reversal = []
for g in markers:
    d = de[g]
    marker_reversal.append({
        "gene": g,
        "logFC": round(d["logFC"], 4),
        "t": round(d["t"], 3),
        "p": d["p"],
        "FDR": d["FDR"],
        "direction_correct": d["logFC"] < 0,   # expected: dupilumab down-regulates type-2 markers
        "significant_FDR05": d["FDR"] < 0.05,
    })

# (c) program score changes under dupilumab (step8 program_summary.csv)
summ_path = os.path.join(BASE, "analysis/step8_expansion/programs/program_summary.csv")
prog = []
with open(summ_path) as f:
    for row in csv.DictReader(f):
        if row["cohort"] == "GSE267331":
            prog.append({
                "program": row["program"],
                "mean_post": float(row["mean_case"]),
                "mean_pre": float(row["mean_ctrl"]),
                "delta_post_minus_pre": round(float(row["diff"]), 6),
                "p_paired_signedrank": float(row["wilcox_p"]),
                "FDR": float(row["FDR"]),
                "n_pairs": int(row["n_case"]),
                "test": row["test"],
            })
prog.sort(key=lambda r: r["p_paired_signedrank"])

# (b) treatment contrasts actually run (from step2 runlog.csv)
runlog = os.path.join(BASE, "analysis/step2_bulk/runlog.csv")
treat_contrasts, other_tx = [], []
with open(runlog) as f:
    for row in csv.DictReader(f):
        c, gse, note = row["contrast"], row["gse"], row["note"]
        if "post_vs_pre" in c and gse == "GSE267331":
            treat_contrasts.append({"gse": gse, "contrast": c,
                "n_case": int(row["n_case"]), "n_control": int(row["n_control"]),
                "method": row["method"], "note": note.strip()})
        elif "stim" in c.lower() or "stimulated" in c:
            other_tx.append({"gse": gse, "contrast": c, "note": note.strip()})

pack = {
    "cohorts": [{
        "accession": "GSE267331",
        "design": "SINUS-52 randomized trial, paired within-subject pre/post nasal brushing; limma paired model",
        "n": {"dupilumab_pairs": 47, "placebo_pairs": 19},
        "contrasts": [c["contrast"] for c in treat_contrasts],
        "contrast_details": treat_contrasts,
        "stats": {
            "dupilumab_DE_FDR05": 4856,
            "placebo_DE_FDR05": 0,
            "module_program_changes": prog,
            "marker_reversal_summary": "7/8 type-2 markers direction-correct (IL33 essentially unchanged, logFC=+0.016, FDR~1.0); 3/8 FDR<0.05 (POSTN, SERPINB2, CLC)",
        },
        "caveats": [
            "Nasal brushing samples, not polyp tissue — cellular context differs from the polyp-resource core.",
            "Readout is marker-level reversal (mRNA), not clinical endpoints (NPS, SNOT-22).",
            "8-gene panel direction check: 7/8 direction-correct per DE file; IL33 essentially flat (+0.016, FDR~1.0).",
            "Placebo arm (19 pairs): 0 genes FDR<0.05 — the reversal is treatment-specific.",
            "Program scores are ssGSEA-style estimates, labeled as estimates, not measured protein."
        ],
    }],
    "marker_reversal": marker_reversal,
    "program_scores_changes": prog,
    "treatment_contrasts_inventory": {
        "clinical_treatment": treat_contrasts,
        "invitro_stimulation_note": other_tx,
        "note": "No steroid, surgery, or other clinical treatment contrasts were run among the 30 bulk contrasts."
    },
}
with open(os.path.join(OUT, "f10_treatment.json"), "w") as f:
    json.dump(pack, f, indent=2)
print("wrote f10_treatment.json")
print("dupilumab contrasts:", len(treat_contrasts), "| invitro:", len(other_tx))
