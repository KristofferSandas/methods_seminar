library(clusterProfiler)
library(edgeR)
library(openxlsx)


organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

# reading in input from deseq2
<<<<<<< Updated upstream
df = read.xlsx("processed_data/DEGs_from_E-MTAB-2523.xlsx", rowNames =TRUE)

# we want the log2 fold change
original_gene_list <- df$logFC
=======
df = read.xlsx("DEGs_from_E-MTAB-2523.xlsx", rowNames =TRUE, colNames = TRUE)



# we want the log2 fold change 
gene_list <- df$logFC
>>>>>>> Stashed changes

# name the vector
names(gene_list) <- rownames(df)



<<<<<<< Updated upstream
# omit any NA values
gene_list<-na.omit(original_gene_list)
=======
# omit any NA values 
#gene_list<-na.omit(original_gene_list)

>>>>>>> Stashed changes

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

#names(gene_list)


# Exctract significant results (padj < 0.05)
#sig_genes_df = subset(df, PValue < 0.05)

# From significant results, we want to filter on log2fold change
#genes <- sig_genes_df$logFC

# Name the vector
names(genes) <- rownames(sig_genes_df)



# omit NA values
#genes <- na.omit(genes)

# filter on min log2fold change (log2FoldChange > 2)
genes <- names(genes)
#genes2 <- stringr::str_to_title(genes) 

#Using bitr function to convert IDs to EntrezIDs

genes_entrez = bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#genes_entrez$ENTREZID
genes_entrez = na.omit(genes_entrez)

#genes_2_entrez = bitr(genes2,fromType = "SYMBOL", toType = "ENTREZID",  OrgDb="org.Hs.eg.db")
#head(genes_2_entrez)

go_enrich2 <- enrichGO(gene = genes,universe = gene_list, OrgDb = organism, keyType = 'SYMBOL',ont = "ALL",readable = T, pvalueCutoff = 0.05, qvalueCutoff = 0.1)
go_enrich

