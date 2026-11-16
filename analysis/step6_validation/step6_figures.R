suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
OUT <- "."
save_fig <- function(p, name, w=8, h=6){
  ggsave(file.path(OUT, paste0(name,".png")), p, width=w, height=h, dpi=300)
  ggsave(file.path(OUT, paste0(name,".pdf")), p, width=w, height=h)
}
## ---- Fig A: marker recovery heatmap ----
mr <- read.csv("marker_recovery.csv", stringsAsFactors=FALSE)
ha <- mr %>% filter(grepl("spatial", layer)==FALSE, layer %in% c("bulk","scrna")) %>%
  mutate(cls = ifelse(layer=="bulk", paste0("Bulk: ", contrast_class), "scRNA: GSE235711 CRSwNP"),
         frac_dir_correct = as.numeric(frac_dir_correct), frac_dir_and_sig = as.numeric(frac_dir_and_sig),
         lbl = ifelse(is.na(frac_dir_correct), "not\nmeasured", sprintf("%.0f%%\n(%.0f%%)", 100*frac_dir_correct, 100*frac_dir_and_sig)),
         fillv = ifelse(is.na(frac_dir_correct), 0.5, frac_dir_correct),
         gene = factor(gene, levels=rev(c("CLC","CCL26","CCL13","IL33","TSLP","POSTN","SERPINB2","ALOX15"))),
         cls = factor(cls, levels=c("Bulk: disease","Bulk: treatment","scRNA: GSE235711 CRSwNP")))
pA <- ggplot(ha, aes(x=cls, y=gene, fill=fillv)) +
  geom_tile(color="grey40") +
  geom_text(aes(label=lbl), size=3.2, color="black") +
  scale_fill_gradient2("direction-correct", low="#d7301f", mid="#fee8c8", high="#1a9850",
                       midpoint=0.5, limits=c(0,1), labels=scales::percent) +
  labs(title="Known-marker recovery across resource layers",
       subtitle="Top %: direction-correct; bottom %: direction-correct AND FDR<0.05\nBulk-disease: 7 contrasts (n>=5/arm, disease-vs-control prior) | treatment: GSE267331 dupilumab | scRNA: GSE235711 CRSwNP",
       x=NULL, y=NULL) + theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=20, hjust=1))
save_fig(pA, "figA_marker_recovery", 10.5, 6.5)
## ---- Fig B: consensus genes ----
cg <- read.csv("consensus_genes.csv", stringsAsFactors=FALSE)
cb <- cg %>% filter(version=="disease_cohorts_10") %>% group_by(direction) %>%
  arrange(desc(n_sig_FDR005), desc(n_consistent)) %>% slice_head(n=15) %>%
  ungroup() %>% arrange(direction, desc(n_sig_FDR005)) %>%
  mutate(gene = factor(gene, levels=unique(gene)))
pB <- ggplot(cb, aes(x=gene, y=frac_consistent, fill=direction)) +
  geom_col() + coord_flip() +
  geom_text(aes(label=paste0(n_sig_FDR005," sig / ",n_consistent,"/",n_cohorts_measured," cohorts")), hjust=1.05, size=2.6) +
  scale_fill_manual(values=c(down="#3182bd", up="#e6550d")) +
  labs(title="Consensus genes (disease-cohort bulk evidence)",
       subtitle="Measured in >=6/10 disease cohorts; direction-consistent in >=80%; FDR<0.05 in >=3\nLabel: sig cohorts / consistent cohorts / measured cohorts",
       x=NULL, y="fraction direction-consistent") + theme_minimal(base_size=11)
save_fig(pB, "figB_consensus", 10, 8.5)
## ---- Fig C: reproduction audit ----
au <- read.csv("reproduction_audit.csv", stringsAsFactors=FALSE)
testable <- au %>% filter(verdict %in% c("reproduced","direction_only","not_reproduced","not_measured")) %>%
  mutate(cohort = paste0(gse, "\nPMID ", pmid),
         verdict = factor(verdict, levels=c("reproduced","direction_only","not_measured","not_reproduced")))
cnt <- testable %>% count(cohort, verdict)
pC <- ggplot(cnt, aes(x=cohort, y=n, fill=verdict)) +
  geom_col() +
  scale_fill_manual("outcome",
    values=c(reproduced="#1a9850", direction_only="#fee08b", not_measured="#d9d9d9", not_reproduced="#d7301f"),
    labels=c("reproduced (dir + FDR<0.05)","direction-only","not measured in outputs","not reproduced")) +
  geom_text(aes(label=n), position=position_stack(vjust=0.5), size=2.8) +
  labs(title="Original-conclusion reproduction audit",
       subtitle="Per-cohort testable claims extracted from the data-generating paper's abstract; tested against our reanalysis tables",
       x=NULL, y="claims") + theme_minimal(base_size=11) +
  theme(axis.text.x=element_text(angle=35, hjust=1))
save_fig(pC, "figC_reproduction", 11, 6)
cat("figures done\n")
