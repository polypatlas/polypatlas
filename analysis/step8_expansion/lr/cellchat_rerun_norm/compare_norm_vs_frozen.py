#!/usr/bin/env python3
"""Compare normalized CellChat rerun vs frozen (unnormalized) results.
All numbers computed from the CSV files on disk. Writes stats JSON for the report."""
import csv, json, math, os
from scipy.stats import spearmanr

ORIG = "analysis/step8_expansion/lr"
NEW = os.path.join(ORIG, "cellchat_rerun_norm")
KEY = ("source", "target", "ligand", "receptor")

def load(path, discol):
    rows = list(csv.DictReader(open(path)))
    out = {}
    for r in rows:
        k = tuple(r[c] for c in KEY)
        out[k] = {"prob_ctrl": float(r["prob_ctrl"]), "prob_dis": float(r[discol]),
                  "delta": float(r["delta"])}
    return out

def ranks(d):
    order = sorted(d, key=lambda k: -abs(d[k]["delta"]))
    return {k: i + 1 for i, k in enumerate(order)}

res = {}
for tag, discol in [("CRSwNP", "prob_dis"), ("CRSsNP", "prob_sn")]:
    old = load(os.path.join(ORIG, f"lr_interactions_{tag}_vs_control.csv"), discol)
    new = load(os.path.join(NEW, f"lr_interactions_{tag}_vs_control_norm.csv"), discol)
    common = set(old) & set(new)
    rho_delta, _ = spearmanr([old[k]["delta"] for k in common], [new[k]["delta"] for k in common])
    rho_absrank, _ = spearmanr([abs(old[k]["delta"]) for k in common], [abs(new[k]["delta"]) for k in common])
    res[tag] = {
        "old_rows": len(old), "new_rows": len(new), "common": len(common),
        "common_frac_of_old": len(common) / len(old), "common_frac_of_new": len(common) / len(new),
        "old_only": len(set(old) - set(new)), "new_only": len(set(new) - set(old)),
        "spearman_delta_common": rho_delta, "spearman_absdelta_common": rho_absrank,
        "sign_agree_common": sum(1 for k in common if (old[k]["delta"] > 0) == (new[k]["delta"] > 0)) / len(common),
    }
    ro, rn = ranks(old), ranks(new)
    # headline pairs named in manuscript §3.7 / SuppFig 1 caption
    headlines = [
        ("tuft", "monocyte", "APP", "CD74"),
        ("tuft", "pDC", "APP", "CD74"),
        ("tuft", "neutrophil", "ANXA1", "FPR1"),
        ("tuft", "pDC", "MIF", "CD74_CXCR4"),
        ("tuft", "monocyte", "ANXA1", "FPR1"),
        ("tuft", "B", "APP", "CD74"),
    ]
    hl = []
    for k in headlines:
        hl.append({"pair": f"{k[0]}->{k[1]} {k[2]}-{k[3]}",
                   "old_delta": old.get(k, {}).get("delta"), "old_rank": ro.get(k),
                   "new_delta": new.get(k, {}).get("delta"), "new_rank": rn.get(k),
                   "in_old": k in old, "in_new": k in new})
    res[tag]["headlines"] = hl
    # top-12 gains / losses old vs new (SuppFig 1 selection)
    def top(d, sign, n=12):
        ks = [k for k in d if (d[k]["delta"] > 0) == (sign > 0)]
        ks.sort(key=lambda k: -d[k]["delta"] if sign > 0 else d[k]["delta"])
        return ks[:n]
    for nm, sign in [("gain", 1), ("loss", -1)]:
        to, tn = top(old, sign), top(new, sign)
        res[tag][f"top12_{nm}_overlap"] = len(set(to) & set(tn))
        res[tag][f"top12_{nm}_old"] = [f"{k[0]}->{k[1]} {k[2]}-{k[3]} d={old[k]['delta']:.3e}" for k in to]
        res[tag][f"top12_{nm}_new"] = [f"{k[0]}->{k[1]} {k[2]}-{k[3]} d={new[k]['delta']:.3e}" for k in tn]
    # largest |delta| rows overall (the frozen table's dominant losses)
    def topabs(d, n=10):
        ks = sorted(d, key=lambda k: -abs(d[k]["delta"]))[:n]
        return [f"{k[0]}->{k[1]} {k[2]}-{k[3]} d={d[k]['delta']:.3e}" for k in ks]
    res[tag]["top10_abs_old"] = topabs(old); res[tag]["top10_abs_new"] = topabs(new)
    # type-2 cytokine ligands presence
    for lig in ["IL33", "TSLP", "IL13", "IL4", "IL5", "CCL26", "CCL13"]:
        res[tag][f"lig_{lig}"] = {"old": sum(1 for k in old if k[2] == lig),
                                  "new": sum(1 for k in new if k[2] == lig)}

print(json.dumps(res, indent=1))
with open(os.path.join(NEW, "comparison_stats.json"), "w") as f:
    json.dump(res, f, indent=1)
