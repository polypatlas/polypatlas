#!/usr/bin/env Rscript
# F6 (v1.6, 2026-11-16): literature support badge from the internal reproduction
# audit v1.7 (85 rows; 53 gene-keyed definitive claims after the round-2
# re-adjudication). Same derivation as f6_lit_v16.R; outputs overwrite
# f6_lit.csv/.js so the website serves the v1.6 badges.
suppressMessages(library(data.table))

aud <- fread("analysis/step6_validation/reproduction_audit_v1.7.csv")
cat("rows:", nrow(aud), "\n")
cat("columns:", paste(colnames(aud), collapse=","), "\n")
cat("verdict levels:", paste(sort(unique(aud$verdict)), collapse=", "), "\n")

aud <- aud[gene != "" & !is.na(gene)]
# Drop non-gene rows: "N/A" (placeholder for a transcriptome-profile claim, R9)
# and every "*_fraction" key (cell-composition claims: R23-mast; R32b
# IDO1pos_fibroblast_fraction and R32c monocyte_fraction; R32a was re-graded
# not_testable in round 2) — none is a gene symbol and this badge file is
# gene-keyed. Documented, not silent. Result: 53 gene-keyed definitive
# claims (22 reproduced / 28 direction_only / 3 not_reproduced) over 47 genes,
# matching denominator_chain_v3.md.
aud <- aud[!gene %in% c("N/A") & !grepl("_fraction$", gene)]
# claims_tested = claims with a definitive verdict (reproduced/direction_only/
# not_reproduced). The audit file also carries not_testable / not_measured /
# skipped rows which are NOT counted.
agg <- aud[, .(
  reproduced = sum(verdict=="reproduced"),
  direction_only = sum(verdict=="direction_only"),
  not_reproduced = sum(verdict=="not_reproduced")
), by=gene]
agg[, claims_tested := reproduced + direction_only + not_reproduced]
agg <- agg[claims_tested > 0]
agg <- agg[, .(gene, claims_tested, reproduced, direction_only, not_reproduced)]
setorder(agg, gene)
fwrite(agg, "analysis/step9_features/f6_lit.csv")
cat("audited genes:", nrow(agg), "\n")
cat("reproduced-any genes:", sum(agg$reproduced>0), "\n")

con <- file("analysis/step9_features/f6_lit.js", "w")
writeLines("var F6LIT = {", con)
n <- nrow(agg)
for (i in seq_len(n)) {
  r <- agg[i]
  writeLines(sprintf('"%s":{t:%d,r:%d,d:%d,n:%d}%s', r$gene, r$claims_tested,
                     r$reproduced, r$direction_only, r$not_reproduced,
                     if (i<n) "," else ""), con)
}
writeLines("};", con); close(con)
cat("wrote f6_lit.csv + f6_lit.js\n")

## ---- round-2 assertions ----------------------------------------------------
stopifnot(sum(agg$claims_tested) == 53, sum(agg$reproduced) == 22,
          sum(agg$direction_only) == 28, sum(agg$not_reproduced) == 3,
          nrow(agg) == 47)
cat("ASSERTIONS PASSED: 53 gene-keyed definitive claims (22/28/3) over 47 genes\n")
