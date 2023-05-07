---
title: "R Notebook"
output: html_notebook
---


```{r echo=T, include=TRUE}
# load required packages
#BiocManager::install("org.Hs.eg.db")
#organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
#library(organism, character.only = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(data.table)

setwd("C:/Users/carlo/OneDrive/Documenti/University/Intership")

# set the path to the folder containing the CSV files
folder_path <- "/Users/carlo/OneDrive/Documenti/University/Intership/RNA/Results/Therapeutic"

# get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = ".csv")

# check if the Ensembl IDs are in a different column name
ensembl_col_name <- "GeneID"
for (DEGs in csv_files) {
  gene_df <- read.csv(file.path(folder_path, DEGs))
if (ensembl_col_name %in% colnames(gene_df)) {
  gene_list <- gene_df$log2FoldChange
  names(gene_list) <- gene_df[[ensembl_col_name]]
  
  # convert Ensembl IDs to KEGG IDs
  ids<- clusterProfiler::bitr(names(gene_list), fromType = "ENSEMBL",
                              toType = "ENTREZID", OrgDb=org.Hs.eg.db)
  
  gene_df2 = gene_df[gene_df[[ensembl_col_name]] %in% ids$ENSEMBL,]
  ids_unique <- unique(ids[!duplicated(ids$ENSEMBL), ])
  gene_df2 <- gene_df[gene_df[[ensembl_col_name]] %in% ids_unique$ENSEMBL, ]
  gene_df2$ENTREZID = ids_unique$ENTREZID
  
  kegg_gene_list <- gene_df2[, c("log2FoldChange", "padj", "ENTREZID")]
  names(kegg_gene_list) <- c("log2FoldChange", "padj", "ENTREZID")
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # Convert kegg_gene_list_fAD to a data table
  write.table(kegg_gene_list, file = paste0(gsub(".csv", "", DEGs), "_kegg.csv"), sep = ",", quote = FALSE, row.names = FALSE)
  print(paste0("KEGG count data has been stored in: ", DEGs, "_kegg.csv"))
  
  kegg_organism = "hsa"
  gene <- ids$ENTREZID 
} else {
  print(paste0("Ensembl IDs not found in file: ", DEGs))
  }
}

```
