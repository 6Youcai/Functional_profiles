#!/usr/bin/env Rscript

library(stringr)
library(DOSE)

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
#
disease <- enrichDO(gene          = gene_enterzid,
                    ont           = "DO",
                    pvalueCutoff  = 1,
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 1,
                    readable      = TRUE)

barplot(disease, showCategory = 20)
dotplot(disease, showCategory = 20)
log_fc <- genelist$logFC
names(log_fc) <- gene_enterzid
cnetplot(disease, 
         categorySize="pvalue",
         foldChange = log_fc,
         fixed = TRUE)
disease <- as.data.frame(disease)
disease$geneID <- str_replace_all(disease$geneID, "/", ";")
save_table(disease, paste0(out_name, ".unfilter.disease.xls"))
#
ncg <- enrichNCG(gene = gene_enterzid,
                 pvalueCutoff = 1,
                 pAdjustMethod = "BH",
                 qvalueCutoff = 1,
                 readable = TRUE)
ncg <- as.data.frame(ncg)[, -1]
ncg$geneID <- str_replace_all(ncg$geneID, "/", ";")
save_table(ncg, paste0(out_name, ".unfilter.cancer.xls"))
# 
# 
# enrichDGN 
# enrichDGNv 
# ChIPseeker 

# upsetplot(disease)
# enrichMap(disease)
