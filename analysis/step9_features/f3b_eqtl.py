#!/usr/bin/env python3
"""F3b step 2: GTEx v8 significant single-tissue eQTL PRESENCE per gene per tissue.
Proxy tissues: Whole_Blood, Lung, Skin_Sun_Exposed_Lower_leg (GTEx has NO nasal).
Batched gencodeId queries; records presence = >=1 significant eQTL row.
Scope genes come from a JSON list file (argv[1]); output appended to checkpoint."""
import json, time, urllib.request, urllib.parse, os, sys
from datetime import datetime, timezone

OUTDIR = os.path.expanduser("analysis/step9_features")
LOG = os.path.join(OUTDIR, "provenance", "f3b_run.log")
TISSUES = ["Whole_Blood", "Lung", "Skin_Sun_Exposed_Lower_leg"]
BATCH = 100

def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    with open(LOG, "a") as f: f.write(line + "\n")
    print(line, flush=True)

symbols = json.load(open(sys.argv[1]))          # list of symbols
resolve = json.load(open(os.path.join(OUTDIR, "provenance", "f3b_resolve_checkpoint.json")))
pairs = [(s, resolve[s]["gencodeId"]) for s in symbols if s in resolve and resolve[s].get("gencodeId")]
log(f"START eqtl_presence n_symbols={len(symbols)} resolved={len(pairs)}")
presence = {s: {t: False for t in TISSUES} for s, _ in pairs}
for tissue in TISSUES:
    for b in range(0, len(pairs), BATCH):
        chunk = pairs[b:b+BATCH]
        ids = ",".join(g for _, g in chunk)
        page, done_rows = 0, 0
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
            done_rows += len(d.get("data", []))
            if page + 1 >= pg.get("numberOfPages", 1):
                break
            page += 1
            time.sleep(0.4)
        nyes = sum(1 for s in presence if presence[s][tissue])
        log(f"{tissue} batch {b//BATCH+1}/{(len(pairs)+BATCH-1)//BATCH}: rows={done_rows} genes_present={nyes}")
        time.sleep(0.4)
out = os.path.join(OUTDIR, "provenance", f"f3b_presence_{os.path.basename(sys.argv[1]).replace('.json','')}.json")
json.dump({"generated_utc": datetime.now(timezone.utc).isoformat(),
           "dataset": "gtex_v8", "tissues": TISSUES,
           "proxy_note": "GTEx has NO nasal tissue; these are proxy tissues only",
           "presence": presence}, open(out, "w"))
log(f"DONE wrote {out}; per-tissue presence: " +
    ", ".join(f"{t}={sum(1 for s in presence if presence[s][t])}" for t in TISSUES))
