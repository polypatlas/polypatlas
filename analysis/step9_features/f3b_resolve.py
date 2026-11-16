#!/usr/bin/env python3
"""F3b step 1: resolve gene symbols -> versioned GENCODE IDs via GTEx v2 reference/gene.
Resolves consensus-UP 825 (always in F3b scope). Checkpointed, resumable."""
import json, time, urllib.request, urllib.parse, os
from datetime import datetime, timezone

OUTDIR = os.path.expanduser("analysis/step9_features")
CKPT = os.path.join(OUTDIR, "provenance", "f3b_resolve_checkpoint.json")
LOG = os.path.join(OUTDIR, "provenance", "f3b_run.log")

def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    with open(LOG, "a") as f: f.write(line + "\n")
    print(line, flush=True)

req_in = json.load(open(os.path.join(OUTDIR, "f2_request.json")))
up = req_in["full_request"]["data"]["upGenes"] if "full_request" in req_in else req_in["upGenes"]
done = {}
if os.path.exists(CKPT):
    done = json.load(open(CKPT))
log(f"START resolve_up n={len(up)} already={len(done)}")
for i, g in enumerate(up):
    if g in done:
        continue
    rec = {"symbol": g, "gencodeId": None, "status": "ok"}
    try:
        url = "https://gtexportal.org/api/v2/reference/gene?geneId=" + urllib.parse.quote(g)
        with urllib.request.urlopen(urllib.request.Request(url), timeout=60) as r:
            d = json.load(r)
        rows = d.get("data", [])
        if rows:
            rec["gencodeId"] = rows[0].get("gencodeId")
        else:
            rec["status"] = "NOT_FOUND"
    except Exception as e:
        rec["status"] = f"ERROR {type(e).__name__} {str(e)[:120]}"
    done[g] = rec
    if (i + 1) % 100 == 0:
        json.dump(done, open(CKPT, "w"))
        log(f"progress {len(done)}/{len(up)}")
    time.sleep(0.4)
json.dump(done, open(CKPT, "w"))
ok = sum(1 for v in done.values() if v["gencodeId"])
log(f"DONE resolved={ok}/{len(up)}")
