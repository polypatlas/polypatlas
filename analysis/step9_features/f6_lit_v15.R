#!/usr/bin/env Rscript
# F6 (v1.5, 2026-11-16): literature support badge from the internal reproduction
# audit v1.5 (85 rows; 56 gene-keyed definitive claims). Same derivation as
# f6_lit.R, which remains the v1.4-era record; outputs overwrite f6_lit.csv/.js
# so the website serves the v1.5 badges.
suppressMessages(library(data.table))

aud <- fread("analysis/step6_validation/reproduction_audit_v1.5.csv")
cat("rows:", nrow(aud), "\n")
cat("columns:", paste(colnames(aud), collapse=","), "\n")
cat("verdict levels:", paste(sort(unique(aud$verdict)), collapse=", "), "\n")

aud <- aud[gene != "" & !is.na(gene)]
# Drop non-gene rows: "N/A" (placeholder for a transcriptome-profile claim, R9)
# and every "*_fraction" key (cell-composition claims: R23-mast in v1.4; R32a
# LY6Dpos_club_fraction, R32b IDO1pos_fibroblast_fraction and R32c
# monocyte_fraction added in v1.5) — none is a gene symbol and this badge file
# is gene-keyed. Documented, not silent. Result: 56 gene-keyed definitive
# claims (25 reproduced / 25 direction_only / 6 not_reproduced) over 49 genes,
# matching denominator_chain_v2.md.
aud <- aud[!gene %in% c("N/A") & !grepl("_fraction$", gene)]
# claims_tested = claims with a definitive verdict (reproduced/direction_only/
# not_reproduced). The audit file also carries 20 non-testable rows
# (not_testable, not_measured, skipped, excluded) which are NOT counted.
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
