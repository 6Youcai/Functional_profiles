#!/usr/bin/env Rscript

library(AnnotationDbi)
library(clusterProfiler)
library(pathview)
library(stringr)

#
to_symbol <- function(ENTREZIDs) {
  enterz <- str_split(ENTREZIDs, "/")[[1]]
  symlobs <- bitr(enterz,
                  fromType = "ENTREZID",
                  toType   = "SYMBOL",
                  OrgDb    = "org.Hs.eg.db")$SYMBOL
  symlobs <- paste(symlobs, collapse = ";")
  return(symlobs)
}

save_table <- function(out_dat, file_name) {
  write.table(out_dat,
              file      = file_name,
              quote     = FALSE,
              sep       = '\t',
              row.names = FALSE,
              col.names = TRUE)
}

##
argv <- commandArgs(TRUE)
DEG <- argv[1]
out_name <- argv[2]
genelist <- read.table(DEG, sep = '\t', header=T)

gene_enterzid <- bitr(genelist$Gene,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = "org.Hs.eg.db")$ENTREZID
# The annotation package, KEGG.db, is not updated since 2012. 
# It’s now pretty old
kk <- enrichKEGG(gene             = gene_enterzid,
                 organism         = "hsa",
                 pvalueCutoff     = 1,
                 pAdjustMethod    = "BH",
                 # not recommended
                 # use_internal_data = TRUE,
                 qvalueCutoff    = 1)

res <- as.data.frame(kk)
res$geneID <- sapply(res$geneID, to_symbol)
save_table(out_dat   = res, 
           file_name = paste0(out_name, ".unfilter.kegg.xls"))

res_filter <- dplyr::filter(res, pvalue <= 0.05, qvalue < 0.2)
save_table(out_dat   = res_filter,
           file_name = paste0(out_name, ".filtered.kegg.xls"))

pdf(paste0(out_name, ".kegg.pdf"), width = 12)
  dotplot(kk, showCategory = 20)
dev.off()

# path view
log_fc <- genelist$logFC
names(log_fc) <- gene_enterzid

dir.create("pathway")
setwd("pathway")
for(i in 1:nrow(res_filter)) {
  pathview(gene.data  = log_fc,
           pathway.id = res_filter$ID[i],
           species    = "hsa",
           # supply your own data 
           # kegg.dir = "kegg/graph",
           limit     = list(gene = max(log_fc), cpd = 1))
}
