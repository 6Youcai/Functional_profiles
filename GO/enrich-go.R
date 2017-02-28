#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(clusterProfiler)
library(AnnotationDbi)

# read database
geneset <- "gene_association.goa_ref_human.gz"
geneset <- read.delim(geneset, header=FALSE, comment.char = "!")
Term2gene <- data.frame(geneset[5], geneset[3])
genesetdesc <- "go-list.txt"
genesetdesc <- read.delim(genesetdesc, header=FALSE)
Term2name <- data.frame(genesetdesc[1], genesetdesc[2])
Term_type <- genesetdesc[, c(1, 3)]
colnames(Term_type) <- c("ID", "Term type")

###################
save_go <- function(ego, deg_number, filtered) {
  res <- as.data.frame(ego)
  res <- merge(res, Term_type)
  res <- dplyr::select(res, ID, geneID, Description, GeneRatio, BgRatio, `Term type`, pvalue, qvalue, Count)
  res$DEG_list <- deg_number
  # go classification
  go_class <- dplyr::select(res, Description, `Term type`, Count, pvalue)
  go_class <- dplyr::arrange(go_class, pvalue)[1:30, ]
  go_class <- dplyr::arrange(go_class, `Term type`)
  go_class$Description <- factor(go_class$Description, levels = go_class$Description)
  pdf(paste0(Out_name, "_", filtered, ".go.classification.xls"))
  ggplot(go_class, aes(x=Description, y=Count, fill = `Term type`)) +
    geom_bar(stat = "identity") +
    labs(x="GO term", y="Number of genes", title="GO classification") +
    coord_flip() +
    theme_bw()
  dev.off()
  # xls
  xls <- res[, -2]
  write.table(xls,
              file = paste0(Out_name, "_", filtered, ".go.xls"),
              quote = FALSE,
              sep = '\t',
              row.names = FALSE,
              col.names = TRUE)
}

draw_DAG <- function(gene_list, filteredd) {
  if(filteredd == "unfilter") {
    pp <- 1
    pp2 <- 1
  } 
  if(filteredd == "filtered") {
    pp <- 0.05
    pp2 <- 0.2
  } 
  gene_enterzid <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")$ENTREZID
  y1 <- enrichGO(gene_enterzid, OrgDb = "org.Hs.eg.db", ont = "MF",
                 pvalueCutoff = pp, pAdjustMethod = "BH",
                 qvalueCutoff = pp2)
  pdf(paste0(Out_name, "_", filteredd, ".go.MF.pdf"))
    try(plotGOgraph(y1))
  dev.off()
  y2 <- enrichGO(gene_enterzid, OrgDb = "org.Hs.eg.db", ont = "BP",
                 pvalueCutoff = pp, pAdjustMethod = "BH",
                 qvalueCutoff = pp2)
  pdf(paste0(Out_name, "_", filteredd, ".go.BP.pdf"))
    try(plotGOgraph(y2))
  dev.off()
  y3 <- enrichGO(gene_enterzid, OrgDb = "org.Hs.eg.db", ont = "CC",
                 pvalueCutoff = pp, pAdjustMethod = "BH",
                 qvalueCutoff = pp2) 
  pdf(paste0(Out_name, "_", filteredd, ".go.CC.pdf"))
    try(plotGOgraph(y3))
  dev.off()
}

####################
argv <- commandArgs(TRUE)
DEG <- argv[1]
Out_name <- argv[2]

gene <- read.table(DEG, header = TRUE)$Gene
gene_num <- length(gene)

# unfilter and filter
x1 <- enricher(gene, TERM2GENE=Term2gene, TERM2NAME=Term2name, pvalueCutoff = 1, qvalueCutoff = 1)
save_go(x1, gene_num, "unfilter")

x2 <- enricher(gene, TERM2GENE=Term2gene, TERM2NAME=Term2name, pvalueCutoff = 0.05, qvalueCutoff = 0.2)
save_go(x2, gene_num, "filtered")
draw_DAG(gene, "filtered")
