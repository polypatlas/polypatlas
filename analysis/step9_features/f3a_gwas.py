#!/usr/bin/env python3
"""F3a: GWAS Catalog v2 — genome-wide significant (p<=5e-8) associations mapped
to scope genes, for in-scope traits (rhinosinusitis / nasal polyps / asthma /
allergic rhinitis / eosinophil). Polite: 1.0 s sleep between requests per worker;
3 workers => ~0.33 s aggregate. Checkpointed per shard; resumable."""
import json, time, urllib.request, urllib.parse, os, sys
from datetime import datetime, timezone

SHARD_FILE = sys.argv[1] if len(sys.argv) > 1 else None
SHARD_TAG = sys.argv[2] if len(sys.argv) > 2 else "main"
FORCE = len(sys.argv) > 3 and sys.argv[3] == "force"
PAGE_SIZE = 100
SLEEP = 1.5

BASE = "https://www.ebi.ac.uk/gwas/rest/api/v2"
KEYWORDS = ["rhinosinusitis", "nasal polyp", "nasal cavity polyp", "asthma", "allergic rhinitis", "eosinophil"]
P_THRESHOLD = 5e-8

OUTDIR = os.path.expanduser("analysis/step9_features")
CKPT = os.path.join(OUTDIR, "provenance", f"f3a_checkpoint_{SHARD_TAG}.json")
LOG = os.path.join(OUTDIR, "provenance", "f3a_run.log")

scope = json.load(open(os.path.join(OUTDIR, "f3_scope_genes_all.json")))
if SHARD_FILE:
    scope = json.load(open(SHARD_FILE))
order_file = os.path.join(OUTDIR, "f3a_gene_order.json")
if not SHARD_FILE and os.path.exists(order_file):
    order = json.load(open(order_file))
    order = [g for g in order if g in set(scope)]
    order += [g for g in scope if g not in set(order)]
    scope = order
done = {}
if not FORCE:
    main_ckpt = os.path.join(OUTDIR, "provenance", "f3a_checkpoint_main.json")
    for ck in ([main_ckpt] if os.path.exists(main_ckpt) else []) + ([CKPT] if os.path.exists(CKPT) and CKPT != main_ckpt else []):
        for d in json.load(open(ck)):
            done[d["gene"]] = d
    # migrate legacy single checkpoint name
    legacy = os.path.join(OUTDIR, "provenance", "f3a_checkpoint.json")
    if os.path.exists(legacy) and not os.path.exists(main_ckpt):
        for d in json.load(open(legacy)):
            done[d["gene"]] = d

def fetch(url, retries=4):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=90) as r:
                return r.status, r.read()
        except Exception as e:
            if attempt == retries - 1:
                raise
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("unreachable")

def trait_ok(a):
    texts = []
    for et in a.get("efo_traits", []) or []:
        texts.append((et.get("efo_trait") or "").lower())
    for rt in a.get("reported_trait", []) or []:
        texts.append(str(rt).lower())
    return any(kw in t for kw in KEYWORDS for t in texts)

def log(msg):
    line = f"{datetime.now(timezone.utc).isoformat(timespec='seconds')} {msg}"
    with open(LOG, "a") as f:
        f.write(line + "\n")
    print(line, flush=True)

todo = [g for g in scope if g not in done]
log(f"START scope={len(scope)} already_done={len(done)} todo={len(todo)}")
sample_saved = 0
for i, gene in enumerate(todo):
    rec = {"gene": gene, "total_mapped": 0, "inscope_n": 0, "traits": [],
           "best_p": None, "best_pmid": None, "pages": 0, "status": "ok"}
    try:
        page, hits = 0, []
        while True:
            q = urllib.parse.urlencode({"mapped_gene": gene, "page": page, "size": PAGE_SIZE})
            status, raw = fetch(f"{BASE}/associations?{q}")
            data = json.loads(raw)
            if sample_saved < 3:
                safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in gene)
                open(os.path.join(OUTDIR, "provenance",
                     f"f3a_sample_gene{safe}_p{page}.json"), "wb").write(raw)
                sample_saved += 1
            pg = data.get("page", {})
            rec["total_mapped"] = pg.get("totalElements", rec["total_mapped"])
            for a in data.get("_embedded", {}).get("associations", []):
                p = a.get("p_value")
                if p is None or p > P_THRESHOLD:
                    continue
                if not trait_ok(a):
                    continue
                traits = sorted({(et.get("efo_trait") or "") for et in (a.get("efo_traits") or []) if et.get("efo_trait")})
                hits.append({"p": p, "traits": traits, "pmid": a.get("pubmed_id")})
            rec["pages"] = page + 1
            if page + 1 >= pg.get("totalPages", 1):
                break
            page += 1
            time.sleep(SLEEP)
        if hits:
            hits.sort(key=lambda h: h["p"])
            rec["inscope_n"] = len(hits)
            rec["traits"] = sorted({t for h in hits for t in h["traits"]})
            rec["best_p"] = hits[0]["p"]
            rec["best_pmid"] = hits[0]["pmid"]
    except Exception as e:
        rec["status"] = f"ERROR {type(e).__name__}: {str(e)[:150]}"
        log(f"{gene}: {rec['status']}")
    done[gene] = rec
    if (i + 1) % 50 == 0 or i + 1 == len(todo):
        json.dump(list(done.values()), open(CKPT, "w"))
        nhits = sum(1 for d in done.values() if d.get("inscope_n"))
        log(f"progress {len(done)}/{len(scope)} genes_queried; with_in_scope_hits={nhits}")
    time.sleep(SLEEP)
json.dump(list(done.values()), open(CKPT, "w"))
log("DONE")
