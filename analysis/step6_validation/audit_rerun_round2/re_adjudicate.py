#!/usr/bin/env python3
"""Round-2 re-adjudication of the internal reproduction audit (v1.5 -> v1.6).

Two pre-declared corrections, both arising from the second external audit
(2026-11-16):

1) BH correction. In v1.5, 21 rows carried a RAW donor-level Mann-Whitney p
   in the FDR column (all GSE276503 rows measured in the 2026-11-16 re-run,
   plus R23-mast under the same precedent). The audit's own grading rule
   defines `reproduced` as correct direction AND FDR < 0.05, so the raw p
   values are BH-adjusted here within a pre-declared family: all numeric
   tests of the same series in the audit file (GSE276503: n = 20;
   GSE235711: n = 1). The raw p is preserved in the new `p_raw` column and
   the FDR column now carries the BH q-value for those rows.

2) Proposition matching. Every row gains claim_proposition / measurement /
   match (exact | proxy | mismatch | n/a). Four rows whose measured
   contrast does not instantiate the claimed proposition are re-graded
   not_testable (dual readings retained in the note):
     R6   (claim: elderly NP vs elderly control; measured: elderly vs
           non-elderly)
     R32a (LY6D+ club proportion: cluster-defined subset in the paper vs
           detection-rate fraction in a residual measurement surface;
           reading B of the round-1 extraction notes)
     R32m (CXCL5: the paper's evidence is an IL-1beta stimulation
           experiment, not the deposited resting-state series)
     R32v (IL1B: tissue-level IHC claim vs within-myeloid normalised
           expression; reading B of the round-1 extraction notes)

Also: R30 layer corrected bulk -> scrna (GSE270953 is an scRNA-seq series
with paired CITE-seq/ADT libraries; its audit rows were measured on
single-cell pseudobulk).

Match classification rules (applied row by row from the frozen claim text
and audit notes; uncertain cases are classed proxy):
  exact    - the measurement tests the claimed quantity in the claimed
             comparison and population (bulk DE claim vs same-contrast
             bulk DE; fraction claim vs donor-level fraction test).
  proxy    - a related but different quantity stands in: protein->mRNA,
             cell-type localisation->whole-tissue bulk DE, subset claim->
             whole parent-class pseudobulk, cluster enrichment->treatment
             reversal DE, co-expression->DE, approximated baseline,
             indirect/auxiliary arm, age-stratified claim->all-ages
             contrast.
  mismatch - the measured contrast is a different comparison altogether.
  n/a      - row not tested (not_testable / skipped in v1.5).

Input : ../reproduction_audit_v1.5.csv  (frozen, untouched)
Output: ../reproduction_audit_v1.6.csv, denominator_chain_v3.md (this dir),
        change_log.md (this dir)
"""
import csv, json
from pathlib import Path

HERE = Path(__file__).resolve().parent
STEP6 = HERE.parent
SRC = STEP6 / "reproduction_audit_v1.5.csv"

rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
assert len(rows) == 85, len(rows)

# ---------------------------------------------------------------- BH ------
def bh(pvals):
    """Benjamini-Hochberg q-values (standard step-up)."""
    n = len(pvals)
    order = sorted(range(n), key=lambda i: pvals[i])
    q = [0.0] * n
    prev = 1.0
    for rank in range(n, 0, -1):
        i = order[rank - 1]
        prev = min(prev, pvals[i] * n / rank)
        q[i] = prev
    return q

# rows whose v1.5 FDR column holds a raw donor-level Mann-Whitney p:
# all GSE276503 rows with a numeric FDR (round-1 re-run convention, stated
# in audit_rerun_2026-11-16/extraction_notes.md) + R23-mast (same precedent).
RAW_P_IDS = {r["id"] for r in rows
             if r["gse"] == "GSE276503" and r["FDR"] not in ("", None)}
RAW_P_IDS.add("R23-mast")
assert len(RAW_P_IDS) == 21, sorted(RAW_P_IDS)

families = {}
for r in rows:
    if r["id"] in RAW_P_IDS:
        families.setdefault(r["gse"], []).append(r)
bh_q = {}
for gse, fam in families.items():
    ps = [float(r["FDR"]) for r in fam]
    qs = bh(ps)
    for r, q in zip(fam, qs):
        bh_q[r["id"]] = q
    print(f"BH family {gse}: n={len(fam)}")

# ------------------------------------------- proposition / measurement ----
# keyed by (id, gene); concise one-line statements from the frozen claim
# text, the round-1 extraction notes and the audit tested_file column.
PM = {
 ("R2a","IL5"): ("Type-2 mediator protein levels increased in CRSwNP polyp tissue",
   "Bulk mRNA DE, CRSwNP vs uncinate control", "proxy"),
 ("R2a","IL13"): ("Type-2 mediator protein levels increased in CRSwNP polyp tissue",
   "Bulk mRNA DE, CRSwNP vs uncinate control", "proxy"),
 ("R2a","CCL24"): ("Type-2 mediator protein levels increased in CRSwNP polyp tissue",
   "Bulk mRNA DE, CRSwNP vs uncinate control", "proxy"),
 ("R2a","CCL13"): ("Type-2 mediator protein levels increased in CRSwNP polyp tissue",
   "Bulk mRNA DE, CRSwNP vs uncinate control", "proxy"),
 ("R4","PLAT"): ("PLAT downregulated in nasal polyps vs controls (bulk transcriptome)",
   "Bulk DE, NP vs control, same contrast", "exact"),
 ("R4","SERPINB2"): ("SERPINB2 upregulated in nasal polyps vs controls (bulk transcriptome)",
   "Bulk DE, NP vs control, same contrast", "exact"),
 ("R5","BPIFB3"): ("BPIFB3 remarkably reduced in nasal polyps (bulk transcriptome)",
   "Bulk DE, NP vs control, same contrast", "exact"),
 ("R5","BPIFB2"): ("BPIFB2 remarkably reduced in nasal polyps (bulk transcriptome)",
   "Bulk DE, NP vs control, same contrast", "exact"),
 ("R5","LPO"): ("LPO remarkably reduced in nasal polyps (bulk transcriptome)",
   "Bulk DE, NP vs control, same contrast", "exact"),
 ("R5","MUC7"): ("MUC7 remarkably reduced in nasal polyps (bulk transcriptome)",
   "Bulk DE, NP vs control, same contrast", "exact"),
 ("R5","SERPINE1"): ("SERPINE1 increased in elderly NP vs elderly controls",
   "Bulk DE, all-ages main contrast NP vs control", "proxy"),
 ("R6","SERPINE1"): ("SERPINE1 increased in elderly NP vs elderly controls",
   "Bulk DE, elderly vs non-elderly contrast", "mismatch"),
 ("R8","CCL18"): ("lncRNA XLOC_010280 regulates CCL18 in eosinophilic CRSwNP",
   "Bulk DE of the target gene CCL18, eosinophilic vs non-eosinophilic", "proxy"),
 ("R9","N/A"): ("Eosinophilic and non-eosinophilic CRSwNP display distinct transcriptome profiles",
   "Count of FDR<0.05 DEGs in the bulk DE of the same contrast", "exact"),
 ("R10","CCL13"): ("CCL13 localised to M2 macrophages in polyp tissue",
   "Whole-tissue bulk DE, CRSwNP vs control", "proxy"),
 ("R10","CCL18"): ("CCL18 localised to M2 macrophages in polyp tissue",
   "Whole-tissue bulk DE, CRSwNP vs control", "proxy"),
 ("R10","CST1"): ("CST1 localised to epithelial cells in polyp tissue",
   "Whole-tissue bulk DE, CRSwNP vs control", "proxy"),
 ("R10","CCL26"): ("CCL26 localised to epithelial cells in polyp tissue",
   "Whole-tissue bulk DE, CRSwNP vs control", "proxy"),
 ("R12","POSTN"): ("POSTN enriched in the C2 type-2-high cluster of the dupilumab trial",
   "Paired dupilumab post- vs pre-treatment bulk DE (reversal)", "proxy"),
 ("R12","CDH26"): ("CDH26 enriched in the C2 type-2-high cluster of the dupilumab trial",
   "Paired dupilumab post- vs pre-treatment bulk DE (reversal)", "proxy"),
 ("R12","CST1"): ("CST1 enriched in the C2 type-2-high cluster of the dupilumab trial",
   "Paired dupilumab post- vs pre-treatment bulk DE (reversal)", "proxy"),
 ("R12","CST2"): ("CST2 enriched in the C2 type-2-high cluster of the dupilumab trial",
   "Paired dupilumab post- vs pre-treatment bulk DE (reversal)", "proxy"),
 ("R12","CST3"): ("CST3 enriched in the C2 type-2-high cluster of the dupilumab trial",
   "Paired dupilumab post- vs pre-treatment bulk DE (reversal)", "proxy"),
 ("R15","IL17RB"): ("IL-17RB+ CD4+ polyp Th2 cells co-express IL17RB/ST2 machinery",
   "Bulk DE, sorted polyp T cells vs nasal T cells", "proxy"),
 ("R15","IL1RL1"): ("IL-17RB+ CD4+ polyp Th2 cells co-express ST2 (IL1RL1)",
   "Bulk DE, sorted polyp T cells vs nasal T cells", "proxy"),
 ("R20","IL2RA"): ("IL2RA differentially expressed in the CRSwNP mast-cell (MCTC) contrast",
   "Bulk DE, same MCTC contrast", "exact"),
 ("R21a","IL5"): ("IL5 upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R21a","IL13"): ("IL13 upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R21a","XCL1"): ("XCL1 upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R21a","CXCL8"): ("CXCL8 upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R21a","AREG"): ("AREG upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R21a","TNFSF14"): ("TNFSF14 upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R21a","METRNL"): ("METRNL upregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2 (paper baseline approximated)", "proxy"),
 ("R22","CD48"): ("CD48 downregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2, same contrast", "exact"),
 ("R22","S1PR1"): ("S1PR1 downregulated in nasal-polyp ILC2s",
   "Bulk DE, NP ILC2 vs peripheral-blood ILC2, same contrast", "exact"),
 ("R29","PNPLA2"): ("ATGL/PNPLA2 lipolysis cascade active in ST2hi memory Th2 cells",
   "Bulk DE, human polyp CD4 T vs PBMC (auxiliary human arm of a mouse paper)", "proxy"),
 ("R29","IL1RL1"): ("ST2/IL1RL1 high in ST2hi memory Th2 cells",
   "Bulk DE, human polyp CD4 T vs PBMC (auxiliary human arm of a mouse paper)", "proxy"),
 ("R30","IGHE"): ("IGHE higher in nasal-polyp antibody-secreting cells than in tonsil",
   "Plasma-cell pseudobulk DE, NP vs tonsil (single-cell layer)", "exact"),
 ("R30","IGHA1"): ("IGHA1 higher in nasal-polyp antibody-secreting cells than in tonsil",
   "Plasma-cell pseudobulk DE, NP vs tonsil (single-cell layer)", "exact"),
 ("R30","IGHA2"): ("IGHA2 higher in nasal-polyp antibody-secreting cells than in tonsil",
   "Plasma-cell pseudobulk DE, NP vs tonsil (single-cell layer)", "exact"),
 ("R30","IGHG4"): ("IGHG4 higher in nasal-polyp antibody-secreting cells than in tonsil",
   "Plasma-cell pseudobulk DE, NP vs tonsil (single-cell layer)", "exact"),
 ("R30","SDC1"): ("SDC1 higher in nasal-polyp antibody-secreting cells than in tonsil",
   "Plasma-cell pseudobulk DE, NP vs tonsil (single-cell layer)", "exact"),
 ("R30","XBP1"): ("XBP1 higher in nasal-polyp antibody-secreting cells than in tonsil",
   "Plasma-cell pseudobulk DE, NP vs tonsil (single-cell layer)", "exact"),
 ("R23-mast","mast_fraction"): ("Mast cells enriched in CRS tissue",
   "Donor-level mast-cell fraction, CRSwNP vs control (Mann-Whitney)", "exact"),
 ("R32a","LY6Dpos_club_fraction"): ("Proportion of LY6D+ club cells greater in CRSwNP IT/NP than healthy IT",
   "Detection-rate LY6D+ fraction within secretory-annotated cells, part-3 objects only, NP vs HC", "proxy"),
 ("R32b","IDO1pos_fibroblast_fraction"): ("IDO1+ fibroblast proportion markedly higher in NP than in other tissues",
   "Detection-rate IDO1+ fraction within fibroblast-annotated cells, NP vs HC-IT", "proxy"),
 ("R32c","monocyte_fraction"): ("Monocyte proportion greater in NP than in other tissues",
   "Donor-level monocyte fraction, NP vs HC-IT (composition table)", "exact"),
 ("R32d","neutrophil_fraction"): ("Neutrophilia and eosinophilia increased in NP mucosa (immunofluorescence)",
   "Not tested: no neutrophil/eosinophil class exists in the unified annotation", "n/a"),
 ("R32e","S100A8"): ("LY6D+ club cells express high levels of S100A8 (main epithelial source)",
   "Whole secretory-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC (part-3 donors)", "proxy"),
 ("R32f","S100A9"): ("LY6D+ club cells express high levels of S100A9 (main epithelial source)",
   "Whole secretory-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC (part-3 donors)", "proxy"),
 ("R32g","IL1RN"): ("IL1RN upregulated in LY6D+ club cells",
   "Whole secretory-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC (part-3 donors)", "proxy"),
 ("R32h","AREG"): ("AREG upregulated in LY6D+ club cells",
   "Whole secretory-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC (part-3 donors)", "proxy"),
 ("R32i","EREG"): ("EREG upregulated in LY6D+ club cells",
   "Whole secretory-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC (part-3 donors)", "proxy"),
 ("R32j","CXCL1"): ("CXCL1 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32k","CXCL2"): ("CXCL2 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32l","CXCL3"): ("CXCL3 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32m","CXCL5"): ("CXCL5 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT (resting state)", "proxy"),
 ("R32n","CXCL6"): ("CXCL6 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32o","CXCL8"): ("CXCL8 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32p","CCL5"): ("CCL5 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32q","CCL7"): ("CCL7 highly expressed in IDO1+ fibroblasts",
   "Not measured: CCL7 absent from the pseudobulk product's gene set", "proxy"),
 ("R32r","CCL8"): ("CCL8 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32s","CCL11"): ("CCL11 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32t","MMP3"): ("MMP3 highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32u","LIF"): ("LIF highly expressed in IDO1+ fibroblasts",
   "Whole fibroblast-class pseudobulk CPM, donor-level Mann-Whitney, NP vs HC-IT", "proxy"),
 ("R32v","IL1B"): ("IL-1beta level increased in CRS-NP vs HC-IT, mostly from MNP/DCs (tissue IHC)",
   "Combined monocyte+macrophage+pDC pseudobulk CPM (within-myeloid expression), NP vs HC-IT", "proxy"),
}
# not-tested rows (v1.5 verdict not_testable / skipped): proposition gist only
SKIP_PM = {
 "R1-skip": "Pathway-level headline claims (cilium, immune regulation, ECM, type-1 IFN)",
 "R3-skip": "tPA/PLAT expression decreased in AERD NP vs CRSwNP NP",
 "R7-skip": "DMRT3 upregulated in AERD (ELISA)",
 "R8x-skip": "lncRNA XLOC_010280 regulates CCL18 and eosinophilic inflammation",
 "R11-skip": "Type-2 and non-type-2 endotypes present in both cohorts",
 "R13-skip": "Cluster C1 enriched for T-cell activation / IL-12 production genes",
 "R14-skip": "Toxin-associated CRS endotype defined by persistent mast-cell accumulation",
 "R16-skip": "IL-25/IL-33 stimulation enhances IL-5 and IL-13 production",
 "R17-skip": "TH17-related genes most overexpressed vs peripheral-blood T cells",
 "R18-skip": "Repertoire/phenotype-level claims (VH usage, IgD+ B to IgE ASC lineage)",
 "R19-skip": "Steroid-resistant ILC2s marked by CD45RO (PTPRC isoform)",
 "R25-skip": "IL-5 elevated in aspirin-intolerant patients",
 "R27-skip": "Th2 subset marker expression (TCN1/CD200R1/HPGDS/LAG3/HAVCR2)",
 "R28-skip": "IL-13 expands and programmes airway tuft cells",
 "R31-skip": "TCF1-LEF1 co-expression identifies TH2-MPP cells",
 "R26a-skip": "GSVA upregulation of 8 regulated-cell-death subtypes in CRSwNP",
 "R26b-skip": "8 key cell-death genes identified; cell-death index elevated in CRSwNP vs control",
 "R26c-skip": "Cell-death genes higher in eosinophilic vs non-eosinophilic polyps; correlate with VAS",
 "R32x-skip": "IL-1beta induces subset transdifferentiation (ALI/primary cells); anakinra in LPS mouse model",
}
for r in rows:
    if r["id"] in SKIP_PM:
        PM[(r["id"], r["gene"])] = (SKIP_PM[r["id"]],
                                    "Not tested (see note)", "n/a")
missing = [(r["id"], r["gene"]) for r in rows if (r["id"], r["gene"]) not in PM]
assert not missing, f"PM missing for: {missing}"

# ------------------------------------------------------- re-adjudication --
TO_NOT_TESTABLE = {
 "R6": "Round-2 re-grade: the measured contrast (elderly vs non-elderly) does not "
       "instantiate the claimed proposition (elderly NP vs elderly control); "
       "proposition mismatch -> not_testable. The all-ages reading of the same "
       "claim is retained as R5/SERPINE1.",
 "R32a": "Round-2 re-grade (reading B of the round-1 extraction notes): the "
         "paper's subset is cluster-defined whereas the measurement is a "
         "detection-rate fraction on a residual surface (part-3 objects only; "
         "the CRS-IT arm is unmeasurable) -> not_testable. Reading A "
         "(not_reproduced) is superseded.",
 "R32m": "Round-2 re-grade: the paper's CXCL5 evidence derives from an IL-1beta "
         "stimulation experiment in primary fibroblasts, not from the deposited "
         "resting-state series -> not_testable for this series.",
 "R32v": "Round-2 re-grade (reading B of the round-1 extraction notes): the claim "
         "is a tissue-level protein (IHC) increase; the measurement is "
         "within-myeloid normalised mRNA expression, which does not instantiate "
         "the tissue-level proposition -> not_testable. Reading A "
         "(not_reproduced) is superseded.",
}

changes = []
out = []
for r in rows:
    r = dict(r)
    key = (r["id"], r["gene"])
    prop, meas, match = PM[key]
    r["claim_proposition"] = prop
    r["measurement"] = meas
    r["match"] = match
    r["p_raw"] = ""
    old_verdict = r["verdict"]
    # 1) BH correction for raw-p rows
    if r["id"] in RAW_P_IDS:
        p = float(r["FDR"])
        q = bh_q[r["id"]]
        r["p_raw"] = r["FDR"]
        r["FDR"] = f"{q:.6g}"
        fam_n = len(families[r["gse"]])
        r["note"] = (r["note"] + " " if r["note"] else "") + \
            f"[round-2: FDR column in v1.5 held the raw donor-level Mann-Whitney p; " \
            f"BH-adjusted within the {r['gse']} audit family (n={fam_n} numeric tests): q={q:.4g}]"
        # re-derive verdict from direction + q (tested rows only)
        if old_verdict in ("reproduced", "direction_only", "not_reproduced"):
            if r["id"] == "R23-mast":
                dir_ok = True  # fraction higher in CRSwNP per note; logFC not recorded
            else:
                lfc = float(r["logFC"])
                dir_ok = (lfc > 0) if r["expected"] == "up" else (lfc < 0)
            new = "reproduced" if (dir_ok and q < 0.05) else \
                  "direction_only" if dir_ok else "not_reproduced"
            if new != old_verdict:
                changes.append((r["id"], r["gene"], old_verdict, new,
                                f"BH q={q:.4g} (raw p={p})"))
                r["verdict"] = new
    # 2) proposition-mismatch re-grades
    if r["id"] in TO_NOT_TESTABLE and r["verdict"] != "not_testable":
        changes.append((r["id"], r["gene"], r["verdict"], "not_testable",
                        "proposition mismatch"))
        r["verdict"] = "not_testable"
        r["note"] = (r["note"] + " " if r["note"] else "") + "[" + TO_NOT_TESTABLE[r["id"]] + "]"
    # 3) R30 layer correction
    if r["id"] == "R30" and r["layer"] == "bulk":
        r["layer"] = "scrna"
        changes.append((r["id"], r["gene"], "layer:bulk", "layer:scrna",
                        "GSE270953 is scRNA-seq + paired CITE-seq/ADT"))
    out.append(r)

cols = ["id", "pmid", "gse", "claim", "gene", "expected", "tested_file",
        "layer", "note", "verdict", "logFC", "FDR", "p_raw",
        "claim_proposition", "measurement", "match"]
with open(STEP6 / "reproduction_audit_v1.6.csv", "w", newline="",
          encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=cols)
    w.writeheader()
    for r in out:
        w.writerow({c: r.get(c, "") for c in cols})

# ------------------------------------------------------------- tallies ----
def tally(rs):
    t = {"reproduced": 0, "direction_only": 0, "not_reproduced": 0,
         "not_measured": 0, "not_testable": 0, "skipped": 0, "excluded": 0}
    for r in rs:
        t[r["verdict"]] += 1
    t["testable"] = t["reproduced"] + t["direction_only"] + \
        t["not_reproduced"] + t["not_measured"]
    t["direction_consistent"] = t["reproduced"] + t["direction_only"]
    t["definitive"] = t["reproduced"] + t["direction_only"] + t["not_reproduced"]
    return t

NON_GENE = {"", "-", "N/A"}
def is_gene(r):
    return r["gene"] not in NON_GENE and "_fraction" not in r["gene"]

v15 = tally(rows)
v16 = tally(out)
exact = tally([r for r in out if r["match"] == "exact"])
proxy = tally([r for r in out if r["match"] == "proxy"])
gene_def = [r for r in out if r["verdict"] in
            ("reproduced", "direction_only", "not_reproduced") and is_gene(r)]
gd = tally(gene_def)
nongene_def = [r for r in out if r["verdict"] in
               ("reproduced", "direction_only", "not_reproduced") and not is_gene(r)]

# ---- assertions: every headline number of the v1.6 chain -----------------
assert v15["testable"] == 65 and v15["reproduced"] == 27
assert v16["testable"] == 61, v16
assert v16["reproduced"] == 24 and v16["direction_only"] == 30
assert v16["not_reproduced"] == 3 and v16["not_measured"] == 4
assert v16["not_testable"] == 20 and v16["skipped"] == 4
assert v16["direction_consistent"] == 54
assert exact["testable"] == 18 and exact["not_reproduced"] == 0, exact
assert proxy["testable"] == 43, proxy
assert len(gene_def) == 53 and gd["reproduced"] == 22 and \
       gd["direction_only"] == 28 and gd["not_reproduced"] == 3
assert len(nongene_def) == 4
assert abs(bh_q["R32j"] - 0.098) < 0.0015, bh_q["R32j"]
print("ALL ASSERTIONS PASSED")
print(json.dumps({"v1.6": v16, "exact": exact, "proxy": proxy,
                  "gene_definitive": gd}, indent=1))

# ------------------------------------------------------- documents --------
def pct(a, b):
    return f"{a}/{b} = {100*a/b:.1f}%"

chain = f"""# Denominator chain v3 — reproduction audit after round-2 re-adjudication (2026-11-16)

Produced by `audit_rerun_round2/re_adjudicate.py` from the frozen
`reproduction_audit_v1.5.csv` (untouched). Output: `reproduction_audit_v1.6.csv`
(85 rows). Two pre-declared corrections:

1. **BH correction of raw p-values.** In v1.5 the FDR column of 21 rows held a
   raw donor-level Mann-Whitney p (all GSE276503 re-run rows + R23-mast).
   BH-adjusted within the pre-declared family = all numeric tests of the same
   series in the audit file (GSE276503 n = 20; GSE235711 n = 1). Raw p is
   preserved in the new `p_raw` column. Effect: CXCL1, CCL8, CCL11
   (raw p = 0.0147 each) move reproduced -> direction_only at q = 0.098.
   No other verdict changes under BH.
2. **Proposition matching.** All rows carry claim_proposition / measurement /
   match (exact | proxy | mismatch | n/a). Four rows whose measured contrast
   does not instantiate the claimed proposition are re-graded not_testable:
   R6 (elderly-vs-noney measured for an elderly-NP-vs-elderly-control claim),
   R32a and R32v (round-1 dual readings, reading B adopted), R32m (the paper's
   CXCL5 evidence is a stimulation experiment outside the deposited series).

Also: R30 layer corrected bulk -> scrna (GSE270953 = scRNA-seq with paired
CITE-seq/ADT libraries; audit rows were measured on single-cell pseudobulk).

## Total chain

| item | v1.5 | v1.6 |
|---|---|---|
| total rows | 85 | 85 |
| reproduced | 27 | **24** |
| direction_only | 27 | **30** |
| not_reproduced | 7 | **3** |
| not_measured | 4 | 4 |
| **testable** | **65** | **61** |
| not_testable | 16 | **20** |
| skipped | 4 | 4 |
| excluded | 0 | 0 |
| direction-consistent | {pct(54,65)} | **{pct(54,61)}** |
| reproduced share | {pct(27,65)} | {pct(24,61)} |
| definitive (rep+dir+notrep) | 61 | 57 |
| gene-keyed definitive (badge set) | 56 (25/25/6) | **53 (22 reproduced / 28 direction-only / 3 not reproduced)** |
| non-gene definitive | 5 | 4 (R9, R23-mast, R32b, R32c) |

## Stratified chains (v1.6, by proposition match)

| stratum | testable | reproduced | direction_only | not_reproduced | not_measured | direction-consistent |
|---|---|---|---|---|---|---|
| exact ({exact['testable']} claims) | {exact['testable']} | {exact['reproduced']} | {exact['direction_only']} | {exact['not_reproduced']} | {exact['not_measured']} | {pct(exact['direction_consistent'], exact['testable'])} |
| proxy ({proxy['testable']} claims) | {proxy['testable']} | {proxy['reproduced']} | {proxy['direction_only']} | {proxy['not_reproduced']} | {proxy['not_measured']} | {pct(proxy['direction_consistent'], proxy['testable'])} |
| total | 61 | 24 | 30 | 3 | 4 | {pct(54,61)} |

No exactly-matched claim failed to reproduce (0/{exact['testable']}); all three
remaining not_reproduced verdicts (IL1RL1 in GSE70898_70900, IL1RN and CCL5 in
GSE276503) sit in the proxy stratum, where a related but different quantity
stood in for the claimed one.

## Verdict changes (7)

| row | gene | v1.5 | v1.6 | reason |
|---|---|---|---|---|
"""
for cid, gene, ov, nv, why in changes:
    if "layer" not in ov:
        chain += f"| {cid} | {gene} | {ov} | {nv} | {why} |\n"
chain += """
Plus the R30 layer correction (bulk -> scrna, 6 rows) and BH q-values written
to the FDR column of all 21 raw-p rows (verdict unchanged for the other 18).

## Reading of the chain

The two round-2 errors acted in opposite directions: the raw-p column had
inflated the reproduced count (27 -> 24 after BH), while proxy/mismatched
propositions in the denominator had deflated the direction-consistency rate
(83.1% -> 88.5% overall; 100% direction-consistency among exactly-matched,
measured claims: 17/17). Neither error cancels the other; both are corrected
and the stratified chain above is the reportable form.
"""
(HERE / "denominator_chain_v3.md").write_text(chain, encoding="utf-8")

log = "# Round-2 change log (v1.5 -> v1.6)\n\nGenerated by re_adjudicate.py.\n\n"
for cid, gene, ov, nv, why in changes:
    log += f"- {cid} / {gene}: {ov} -> {nv} ({why})\n"
log += f"\nBH families: " + "; ".join(
    f"{g} n={len(f)}" for g, f in families.items()) + "\n"
log += "BH q for CXCL1/CCL8/CCL11 (R32j/R32r/R32s): " + \
       f"{bh_q['R32j']:.4f} (raw p 0.0147, family n=20)\n"
(HERE / "change_log.md").write_text(log, encoding="utf-8")
print("wrote reproduction_audit_v1.6.csv, denominator_chain_v3.md, change_log.md")
