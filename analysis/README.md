# PolypAtlas analysis code

This directory contains the analysis code that produced every result in the
PolypAtlas resource and manuscript (release v1.8). It is released under the
MIT licence (see `../LICENSE`). The website build code is in `../build/`.

## Conventions

- Run all scripts **from the project root** (the directory that contains
  `analysis/`); paths inside the scripts are project-root-relative.
- Inputs are public data: GEO series listed in `../site/data/curation_table.csv`
  (downloaded from NCBI GEO), PubMed abstracts/full texts for the reproduction
  audit, and the GWAS Catalog / GTEx Portal APIs for the annotation layer.
- Intermediate and final result files are not duplicated here; they are in the
  Zenodo release snapshot and in `../site/data/`.

## Pipeline order

1. **Census and bulk reanalysis** — `step2_bulk/`
   `01*_fetch_metadata.R` (GEO metadata), `02_prepare.R`, `02_de.R` /
   `03_de.R` (per-cohort differential expression; DESeq2 for counts, limma
   for arrays/normalised data; probes collapsed to genes by minimum p-value),
   `03/04_assemble_evidence.R` (evidence matrices), `04/05_controls.R`
   (positive/negative controls), `rebuild_runlog.R`.
   - `sensitivity_probe_collapse/` — post-hoc probe-collapse sensitivity
     (maximum-variance rule): `sensitivity_de.R`, `consensus_sensitivity.py`
     (also verifies the consensus rule against the frozen consensus table).
   - `did_dupilimab/` — exploratory difference-in-differences
     (treatment × time) analysis of the dupilumab cohort: `did_dupilimab.R`.
   - `sensitivity_strict_tissue/` — strict whole-tissue sensitivity:
     `strict_consensus.py` (consensus rule on the 7 whole-tissue cohorts;
     primary-contrast inventory) and `strict_meta.R` (strict-set
     meta-analysis).
2. **Single-cell processing** — `step3_scrna/`
   `dl.py` (download), `pipeline.py` / `proc_*.py` (per-dataset processing and
   unified annotation), `stage2*.py`, `fin_276503.py`, `pb_de.R`
   (donor-level pseudobulk differential expression, edgeR).
   - `rebuild_composition_v18.py` — release-v1.8 rebuild of
     `composition_all.csv` from the annotated objects: restores the
     `secretory_club` label that an earlier assembly had collapsed into
     `secretory` (asserted one row per donor and class, fractions summing
     to 1, and exact agreement with the previous table after summing the
     two classes back); `fig5_stats_v18.R` — the Figure 5 composition
     statistics (Kruskal–Wallis, Benjamini–Hochberg FDR across the 19
     classes, Hodges–Lehmann effect estimates).
3. **Spatial layer** — `step4_spatial/pipeline.R` (+ `fix1.R`).
4. **Validation and the reproduction audit** — `step6_validation/`
   `step6_figures.R` (validation figures from the frozen result tables).
   The audit itself was assembled from verbatim claim extraction against
   PubMed records; the two adjudication rounds that produced the current
   audit file are fully scripted:
   - `audit_rerun_2026-11-16/` — re-run for the two series whose deposited
     publication links had been wrong (`measure_rerun.py`,
     `measure_rerun_v2.py`);
   - `audit_rerun_round2/re_adjudicate.py` — round-2 re-adjudication:
     Benjamini–Hochberg correction of the raw Mann–Whitney p-values that the
     first re-run had stored in the FDR column, proposition-match
     stratification (exact / proxy / mismatch), and the v1.6 denominator
     chain. Input: `reproduction_audit_v1.5.csv`; output:
     `reproduction_audit_v1.6.csv`.
5. **Analysis layer** — `step8_expansion/`
   `enrichment/run_enrichment.R` (ORA + GSEA), `programs/run_programs.R`
   (ssGSEA program scores), `deconv/run_deconv.R` (MCP-counter-style
   geometric-mean abundance estimates), `lr/` (CellChat: `cellchat_run.R`,
   `cellchat_snp.R`, `cellchat_compare.R`; the library-size-normalised
   re-run that supersedes them in v1.6 is in `lr/cellchat_rerun_norm/`,
   with its comparison report).
6. **Added-value features** — `step9_features/`
   `f1_meta.R` (cross-cohort DerSimonian–Laird meta-analysis),
   `meta_v2/meta_v2.R` (τ², prediction intervals, leave-one-cohort-out;
   asserts equality with the frozen `f1_meta.csv`), `f3a_gwas.py` /
   `f3b_*.py` / `f3_finalize.py` (genetic annotation), `f4_coexpr_*.R`
   (co-expression), `f6_lit*.R` (literature-support badges from the audit;
   `f6_lit_v16.R` is the current version), `f10_pack.py` / `f10_figure.R`
   (treatment-response page data and figure).
7. **Manuscript figures** — `figures/`
   `make_fig*_v3.R` and `make_suppfig1_v3.R`, drawn under the project's
   figure style v3 (`figstyle_v3.R`, included). Numbers in every figure are
   read from the result files; several scripts assert their headline totals.

## Environment

- R 4.5.3 with (versions used for v1.6): DESeq2 1.50.2, limma 3.66.0,
  edgeR 4.8.2, clusterProfiler 4.18.4, GSVA 2.4.9, CellChat 1.6.1,
  data.table 1.18.6.1, ggplot2 4.0.3.
- Python 3.12 with scanpy 1.12.4, anndata 0.13.4, numpy 2.5.4,
  pandas 3.0.6 (single-cell processing).
- The website is built with `../build/build.py` (Python 3, no dependencies
  beyond the standard library).
