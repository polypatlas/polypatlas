#!/usr/bin/env python3
"""F3b extras: for F3a-hit genes NOT in the consensus-UP set (i.e. DOWN or
audit-only genes with GWAS in-scope hits), resolve GENCODE IDs and query
GTEx v8 eQTL presence in the 3 proxy tissues. Run AFTER F3a completes."""
import json, time, urllib.request, urllib.parse, os, sys
from datetime import datetime, timezone

OUTDIR = os.path.expanduser("analysis/step9_features")
PROV = os.path.join(OUTDIR, "provenance")
LOG = os.path.join(PROV, "f3b_run.log")
TISSUES = ["Whole_Blood", "Lung", "Skin_Sun_Exposed_Lower_leg"]

def log(m):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {m}"
    with open(LOG, "a") as f: f.write(line + "\n")
    print(line, flush=True)

# F3a hits (merge ALL shard checkpoints)
gwas = {}
for fn in sorted(os.listdir(PROV)):
    if fn.startswith("f3a_checkpoint_") and fn.endswith(".json"):
        for d in json.load(open(os.path.join(PROV, fn))):
            gwas[d["gene"]] = d
hit_genes = {g for g, d in gwas.items() if d.get("inscope_n")}
req = json.load(open(os.path.join(OUTDIR, "f2_request.json")))
up = set(req["full_request"]["data"]["upGenes"])
extra = sorted(hit_genes - up)
log(f"extras: {len(extra)} F3a-hit genes outside consensus-UP")
resolve = json.load(open(os.path.join(PROV, "f3b_resolve_checkpoint.json")))
pairs = []
for s in extra:
    if s in resolve and resolve[s].get("gencodeId"):
        pairs.append((s, resolve[s]["gencodeId"]))
        continue
    try:
        url = "https://gtexportal.org/api/v2/reference/gene?geneId=" + urllib.parse.quote(s)
        with urllib.request.urlopen(urllib.request.Request(url), timeout=60) as r:
            rows = json.load(r).get("data", [])
        gid = rows[0]["gencodeId"] if rows else None
        resolve[s] = {"symbol": s, "gencodeId": gid,
                      "status": "ok" if gid else "NOT_FOUND"}
        if gid:
            pairs.append((s, gid))
        time.sleep(0.4)
    except Exception as e:
        resolve[s] = {"symbol": s, "gencodeId": None, "status": f"ERROR {e}"[:120]}
json.dump(resolve, open(os.path.join(PROV, "f3b_resolve_checkpoint.json"), "w"))
log(f"resolved {len(pairs)}/{len(extra)} extras")
presence = {s: {t: False for t in TISSUES} for s, _ in pairs}
for tissue in TISSUES:
    for b in range(0, len(pairs), 100):
        chunk = pairs[b:b + 100]
        page = 0
        while True:
            q = urllib.parse.urlencode([("gencodeId", g) for _, g in chunk] +
                                       [("tissueSiteDetailId", tissue),
                                        ("datasetId", "gtex_v8"),
                                        ("page", page), ("itemsPerPage", 10000)])
            url = "https://gtexportal.org/api/v2/association/singleTissueEqtl?" + q
            with urllib.request.urlopen(urllib.request.Request(url), timeout=120) as r:
                d = json.load(r)
            for row in d.get("data", []):
                sym = (row.get("geneSymbolUpper") or row.get("geneSymbol") or "").upper()
                if sym in presence and row.get("tissueSiteDetailId") == tissue:
                    presence[sym][tissue] = True
            pg = d.get("paging_info", {})
            if page + 1 >= pg.get("numberOfPages", 1):
                break
            page += 1
            time.sleep(0.4)
        time.sleep(0.4)
    log(f"{tissue}: present={sum(1 for s in presence if presence[s][tissue])}")
json.dump({"generated_utc": datetime.now(timezone.utc).isoformat(),
           "dataset": "gtex_v8", "tissues": TISSUES,
           "proxy_note": "GTEx has NO nasal tissue; proxy tissues only",
           "presence": presence},
          open(os.path.join(PROV, "f3b_presence_extras.json"), "w"))
log("extras DONE")
