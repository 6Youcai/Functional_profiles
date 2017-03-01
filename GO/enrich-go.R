#!/usr/bin/env Rscript

library(dplyr)
library(stringr)
library(clusterProfiler)
library(AnnotationDbi)
library(ggplot2)

# term and gene
geneset <- "data/goa_human.gaf.gz"
geneset <- read.delim(geneset, 
                      header = FALSE, 
                      comment.char = "!")
Term2gene <- dplyr::select(geneset, V5, V3) # 
# term and name
genesetdesc <- "data/go_namespace.txt"
genesetdesc <- read.delim(genesetdesc, 
                          header = FALSE, 
                          sep = '\t')
Term2name <- dplyr::select(genesetdesc, V1, V2) # 
# term and type
Term2type <- dplyr::select(genesetdesc, V1, V3) #
colnames(Term2type) <- c("ID", "Term type")

# read argv
argv <- commandArgs(TRUE)
DEG <- argv[1]
Out_name <- argv[2]
gene <- read.table(DEG, header = TRUE)$Gene
gene_num <- length(gene)

###
egmt <- enricher(gene, 
                 TERM2GENE=Term2gene, 
                 TERM2NAME=Term2name, 
                 pvalueCutoff = 1, 
                 qvalueCutoff = 1)
res <- as.data.frame(egmt)
res <- merge(res, Term2type, by = 'ID')
res <- dplyr::select(res, -p.adjust)
res$DEG_list <- gene_num

# go classification
go_class <- dplyr::select(res, Description, `Term type`, Count, pvalue)
if(nrow(go_class)>= 30) {
  go_class <- dplyr::arrange(go_class, pvalue)[1:30, ]
}
go_class <- dplyr::arrange(go_class, `Term type`)
go_class$Description <- factor(go_class$Description, 
                               levels = go_class$Description)
pdf(paste0(Out_name, ".go.classification.pdf"),
    height = 12, width = 10)
  ggplot(go_class, 
         aes(x=Description, y=Count, fill = `Term type`)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(x="", y="Number of genes", title="GO classification") +
    theme(axis.text.x=element_text(angle=45, size = 12, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "bottom")
dev.off()
# xls
xls <- dplyr::select(res, -geneID)
write.table(xls,
            file = paste0(Out_name, "_unfilter", ".go.xls"),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = TRUE)
# xls_filter
xls_filter <- filter(xls, pvalue <= 0.05, qvalue < 0.2)
write.table(xls,
            file = paste0(Out_name, "_filtered", ".go.xls"),
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = TRUE)

# draw_DAG
gene_enterzid <- bitr(gene, 
                      fromType="SYMBOL", 
                      toType="ENTREZID", 
                      OrgDb="org.Hs.eg.db")$ENTREZID

y1 <- enrichGO(gene_enterzid,
               ### universe = 
               OrgDb = "org.Hs.eg.db", 
               ont = "MF",
               readable = TRUE,
               pvalueCutoff = 0.05, 
               pAdjustMethod = "BH",
               qvalueCutoff = 0.2)
y2 <- enrichGO(gene_enterzid, 
               OrgDb = "org.Hs.eg.db", 
               ont = "BP",
               readable = TRUE,
               pvalueCutoff = 0.05, 
               pAdjustMethod = "BH",
               qvalueCutoff = 0.2)
y3 <- enrichGO(gene_enterzid, 
               OrgDb = "org.Hs.eg.db", 
               ont = "CC",
               readable = TRUE,
               pvalueCutoff = 0.05, 
               pAdjustMethod = "BH",
               qvalueCutoff = 0.2) 

pdf(paste0(Out_name, "_filteredd.go.MF.pdf"))
  plotGOgraph(y1)
dev.off()

pdf(paste0(Out_name, "_filteredd.go.BP.pdf"))
  plotGOgraph(y2)
dev.off()

pdf(paste0(Out_name, "_filteredd.go.CC.pdf"))
  plotGOgraph(y3)
dev.off()
