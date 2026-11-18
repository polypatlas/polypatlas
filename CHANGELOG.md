# PolypAtlas changelog

Version history of the resource releases. Each release is archived as a fixed snapshot on Zenodo (all-versions DOI: 10.5281/zenodo.23448779); archived snapshots are never edited silently. The manuscript's Data availability section points here and to the version notes of each Zenodo release.

## v1.8 (18 November 2026) — Zenodo DOI 10.5281/zenodo.24329815

Release v1.8 (18 November 2026) supersedes v1.7: it restores the secretory_club label in the single-cell composition table (merged into 'secretory' in earlier releases, which had double-counted donors in the Figure 5 test), recomputes the Figure 5 statistics with FDR control and effect sizes (Supplementary Data 13), corrects the GSE156285 publication chain in the curation table and Table 1 (the data-generating paper [55] restored; the 2025 secondary analysis [50] recorded in its actual role), removes an unsupported 'paratumoral' description of the GSE194282 controls (uncinate tissue from normal healthy controls, per the source paper), renumbers the marker-panel supplement (Table 6 → Table 5) and moves the audit re-adjudication history into Supplementary Note 1;

*Text revision after the v1.8 release (18 November 2026; manuscript package v12, resource unchanged): the reproduction-audit wording was aligned with the §2.6 grading definitions (16/16 measured exactly-matched claims direction-consistent, 9/16 fully reproduced, none direction-inconsistent); the consensus was renamed the "primary disease-contrast consensus" in the manuscript and in the download label above (archived file names unchanged); the negative-results section now names both zero-DEG primary contrasts (GSE69093 and GSE36830); two publication-link errors were re-attributed to our own curation/extraction files; the §2.6 badge denominator chain was completed (56 definitive verdicts minus three gene-less cell-fraction claims = 53 gene-keyed claims); the abstract was shortened to 300 words; and this changelog was created, receiving the release history formerly carried in the manuscript's Data availability section.*

## v1.7 (17 November 2026) — Zenodo DOI 10.5281/zenodo.24300458

Release v1.7 (17 November 2026) supersedes v1.6: it re-adjudicates one study-level audit claim out of the gene-level chain (audit v1.7; 60 testable claims, §2.6), corrects the GSE198950 tissue field in the curation table and Table 1, restores five missing Table 1 citation brackets and adds reference [54], and adds the GSE198950-exclusion and GSE36830 donor-handling sensitivity analyses (Supplementary Data 11 and 12); no primary result changed in either round except as itemised here and in the release changelog.

## v1.6 (16 November 2026) — Zenodo DOI 10.5281/zenodo.24283166

v1.6 (16 November 2026) corrected the GSE270953 design classification (scRNA-seq with paired CITE-seq/ADT libraries; no bulk arm), replaced the reproduction audit with the round-2 re-adjudicated file (Benjamini–Hochberg-corrected single-cell verdicts, proposition-match stratification; 61 testable claims), added the strict whole-tissue consensus and meta-analysis sensitivity tables and the library-size-normalised CellChat re-run, and shipped the analysis code with the snapshot (see Code availability).
