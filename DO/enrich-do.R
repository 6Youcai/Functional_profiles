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

genelist <- read.table(DEG, sep = '\t', header=T)$Gene
gene_enterzid <- bitr(genelist,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = "org.Hs.eg.db")$ENTREZID
disease <- enrichDO(gene          = gene_enterzid,
                    ont           = "DO",
                    pvalueCutoff  = 1,
                    pAdjustMethod = "BH",
                    qvalueCutoff  = 1,
                    readable      = TRUE)
disease <- as.data.frame(disease)
disease$geneID <- str_replace_all(disease$geneID, "/", ";")
save_table(out_dat   = disease, file_name = paste0(out_name, ".unfilter.disease.xls"))
