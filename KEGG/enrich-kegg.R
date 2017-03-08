#!/usr/bin/env Rscript

library(AnnotationDbi)
library(clusterProfiler)
library(pathview)
library(stringr)
library(dplyr)

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

symbol_enterzid <- bitr(genelist$Gene,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = "org.Hs.eg.db")
# The annotation package, KEGG.db, is not updated since 2012. 
# It’s now pretty old
kk <- enrichKEGG(gene             = symbol_enterzid$ENTREZID,
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

# path vis
# while bitr some gene symbol will be faild
# so use merge to ensure the ENTREZID and logFC 
# one by one correctly
colnames(symbol_enterzid) <- c("Gene", "ENTREZID")
fc <- merge(symbol_enterzid, genelist, by = "Gene")
fc_fixed <- fc$logFC
names(fc_fixed) <- fc$ENTREZID

dir.create("pathway")
setwd("pathway")
for(i in 1:nrow(res_filter)) {
  pathview(gene.data  = fc_fixed,
           pathway.id = res_filter$ID[i],
           species    = "hsa",
           # supply your own data
           # kegg.dir = "kegg/graph",
           limit     = list(gene = max(fc_fixed), cpd = 1))
}
