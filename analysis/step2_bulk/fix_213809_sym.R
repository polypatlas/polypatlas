suppressPackageStartupMessages({ library(org.Hs.eg.db) })
OUT <- "analysis/step2_bulk"
sp <- readRDS(file.path(OUT, "prep", "GSE213809.rds"))
ids <- sub("\\.[0-9]+$", "", rownames(sp$expr))
sym <- suppressWarnings(AnnotationDbi::mapIds(org.Hs.eg.db, keys = ids, column = "SYMBOL",
                                              keytype = "ENSEMBL", multiVals = "first"))
names(sym) <- rownames(sp$expr)
cat("mapped:", sum(!is.na(sym)), "/", length(sym), "\n")
sp$symbols <- sym
saveRDS(sp, file.path(OUT, "prep", "GSE213809.rds"))
cat("prep updated\n")
