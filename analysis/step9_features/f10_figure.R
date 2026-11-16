#!/usr/bin/env Rscript
# F10 figure: paired pre/post program scores, GSE267331 dupilumab arm (47 pairs).
# Real per-sample ssGSEA-style program scores from step8; pairing from prep meta.
suppressMessages(library(data.table))

ps <- fread("analysis/step8_expansion/programs/program_scores.csv")
x <- readRDS("analysis/step2_bulk/prep/GSE267331.rds")
meta <- as.data.table(x$meta)

d <- ps[cohort=="GSE267331"]
d <- merge(d, meta[, .(sample, subject, time, treatment)], by="sample")
d <- d[treatment=="Dupilumab" & time %in% c("Pre-Treatment","Post-Treatment")]
d[, time := factor(time, levels=c("Pre-Treatment","Post-Treatment"))]
d <- dcast(d, subject + program ~ time, value.var="score")
setnames(d, c("Pre-Treatment","Post-Treatment"), c("pre","post"))
d <- d[!is.na(pre) & !is.na(post)]
cat("paired subjects:", length(unique(d$subject)), "\n")

summ <- fread("analysis/step8_expansion/programs/program_summary.csv")
summ <- summ[cohort=="GSE267331"]
progs <- summ[FDR < 0.05][order(wilcox_p)]$program
cat("FDR<0.05 programs:", paste(progs, collapse=", "), "\n")

labs <- sapply(progs, function(p) {
  r <- summ[program==p]
  sprintf("%s\npost-pre=%.3f, p=%.1e", p, r$diff, r$wilcox_p)
})

png("analysis/step9_features/f10_dupilumab_paired_scores.png",
    width=10, height=7, units="in", res=300)
par(mfrow=c(2,3), mar=c(4.2,4.5,3.2,1), mgp=c(2.6,0.8,0))
for (i in seq_along(progs)) {
  s <- d[program==progs[i]]
  plot(1:2, range(c(s$pre, s$post)), type="n", xaxt="n",
       xlab="", ylab="program score", main=labs[i], cex.main=0.95,
       xlim=c(0.6, 2.4))
  axis(1, at=1:2, labels=c("Pre","Post"))
  for (k in seq_len(nrow(s))) {
    col <- if (s$post[k] < s$pre[k]) "#2E7D32AA" else "#C62828AA"
    lines(1:2, c(s$pre[k], s$post[k]), col=col, lwd=0.7)
  }
  points(rep(1, nrow(s)), s$pre, pch=19, cex=0.6, col="#1B5E20")
  points(rep(2, nrow(s)), s$post, pch=19, cex=0.6, col="#0D47A1")
  points(1:2, c(mean(s$pre), mean(s$post)), pch=18, cex=1.6, col="black")
  lines(1:2, c(mean(s$pre), mean(s$post)), lwd=2.2, col="black")
}

pdf("analysis/step9_features/f10_dupilumab_paired_scores.pdf", width=10, height=7)
par(mfrow=c(2,3), mar=c(4.2,4.5,3.2,1), mgp=c(2.6,0.8,0))
for (i in seq_along(progs)) {
  s <- d[program==progs[i]]
  plot(1:2, range(c(s$pre, s$post)), type="n", xaxt="n",
       xlab="", ylab="program score", main=labs[i], cex.main=0.95,
       xlim=c(0.6, 2.4))
  axis(1, at=1:2, labels=c("Pre","Post"))
  for (k in seq_len(nrow(s))) {
    col <- if (s$post[k] < s$pre[k]) "#2E7D32AA" else "#C62828AA"
    lines(1:2, c(s$pre[k], s$post[k]), col=col, lwd=0.7)
  }
  points(rep(1, nrow(s)), s$pre, pch=19, cex=0.6, col="#1B5E20")
  points(rep(2, nrow(s)), s$post, pch=19, cex=0.6, col="#0D47A1")
  points(1:2, c(mean(s$pre), mean(s$post)), pch=18, cex=1.6, col="black")
  lines(1:2, c(mean(s$pre), mean(s$post)), lwd=2.2, col="black")
}
dev.off()
cat("wrote f10_dupilumab_paired_scores.png + .pdf\n")
