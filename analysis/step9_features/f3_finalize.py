#!/usr/bin/env python3
"""F3 finalize: build f3_genetic.csv + f3_genetic.js from F3a checkpoint and
F3b presence files. Only genes with >=1 hit (GWAS in-scope association or
GTEx eQTL presence in a proxy tissue) are included."""
import json, csv, os
from datetime import datetime, timezone

OUTDIR = os.path.expanduser("analysis/step9_features")
PROV = os.path.join(OUTDIR, "provenance")
PROXY = "GTEx has NO nasal tissue; proxy tissue only"

gwas = {}
for fn in sorted(os.listdir(PROV)):
    if fn.startswith("f3a_checkpoint_") and fn.endswith(".json"):
        for d in json.load(open(os.path.join(PROV, fn))):
            gwas[d["gene"]] = d
gwas = list(gwas.values())
presence_files = [os.path.join(PROV, f) for f in os.listdir(PROV)
                  if f.startswith("f3b_presence_")]
pres = {}
for f in presence_files:
    d = json.load(open(f))
    for sym, tmap in d["presence"].items():
        pres.setdefault(sym.upper(), {})
        pres[sym.upper()].update(tmap)

KEYWORDS = ["rhinosinusitis", "nasal polyp", "nasal cavity polyp", "asthma",
            "allergic rhinitis", "eosinophil"]

def inscope_traits(traits):
    return [t for t in traits
            if any(kw in t.lower() for kw in KEYWORDS)]

rows = []
for x in gwas:
    if x.get("inscope_n"):
        it = inscope_traits(x["traits"])
        rows.append({"gene": x["gene"], "source": "GwasCatalog",
                     "trait_or_tissue": "; ".join(it) if it else "; ".join(x["traits"]),
                     "n_associations": x["inscope_n"],
                     "best_p": x["best_p"], "study_pmid": x["best_pmid"],
                     "proxy_note": ""})
for sym, tmap in sorted(pres.items()):
    for t, on in tmap.items():
        if on:
            rows.append({"gene": sym, "source": "GTEx",
                         "trait_or_tissue": t.replace("_", " ") + " (proxy)",
                         "n_associations": 1, "best_p": "", "study_pmid": "",
                         "proxy_note": PROXY + "; presence = >=1 significant GTEx v8 single-tissue eQTL"})

with open(os.path.join(OUTDIR, "f3_genetic.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["gene", "source", "trait_or_tissue",
                                      "n_associations", "best_p", "study_pmid", "proxy_note"])
    w.writeheader(); w.writerows(rows)

js = {}
for r in rows:
    g = r["gene"]
    js.setdefault(g, []).append({"source": r["source"],
                                 "trait": r["trait_or_tissue"],
                                 "bestp": r["best_p"]})
with open(os.path.join(OUTDIR, "f3_genetic.js"), "w") as f:
    f.write("// PolypAtlas feature F3: genetic support badges.\n")
    f.write(f"// GWAS Catalog v2 (p<=5e-8, in-scope traits) + GTEx v8 eQTL presence (PROXY tissues: GTEx has NO nasal tissue).\n")
    f.write(f"// Generated {datetime.now(timezone.utc).isoformat(timespec='seconds')} UTC. Only genes with >=1 hit.\n")
    f.write("var F3GENETIC = " + json.dumps(js, ensure_ascii=False) + ";\n")

n_gwas_genes = sum(1 for x in gwas if x.get("inscope_n"))
n_gtex_genes = sum(1 for s, t in pres.items() if any(t.values()))
print(f"rows={len(rows)} gwas_genes={n_gwas_genes} gtex_genes_with_presence={n_gtex_genes} js_genes={len(js)}")
print(f"gwas coverage: {len(gwas)} identifiers queried")
