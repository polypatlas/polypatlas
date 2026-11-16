#!/usr/bin/env python3
"""Download GEO series tars via the geo/download endpoint (individual files 404 for most series)."""
import os, sys, subprocess, re
BASE = os.path.expanduser("analysis/step3_scrna/raw")

def dl(acc, expect_mb):
    d = os.path.join(BASE, acc); os.makedirs(d, exist_ok=True)
    out = os.path.join(d, f"{acc}_RAW.tar")
    if os.path.exists(out) and os.path.getsize(out) > 1024:
        print(f"{acc}: already present ({os.path.getsize(out)//1048576} MB)"); return True
    url = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={acc}&format=file"
    print(f"{acc}: downloading (~{expect_mb} MB)...", flush=True)
    r = subprocess.run(["curl","-sL","-o",out,"-w","%{http_code} %{size_download}",url,"--max-time","3600"],
                       capture_output=True, text=True)
    code, size = r.stdout.strip().split()
    ok = code == "200" and int(size) > 1024
    print(f"{acc}: http={code} bytes={size} -> {'OK' if ok else 'FAIL'}", flush=True)
    return ok

if __name__ == "__main__":
    jobs = [("GSE324445",111),("GSE196169",544),("GSE156285",350),("GSE249159",420),
            ("GSE255544",118),("GSE235711",580),("GSE276503",1282)]
    only = sys.argv[1:]
    for acc, mb in jobs:
        if only and acc not in only: continue
        dl(acc, mb)
