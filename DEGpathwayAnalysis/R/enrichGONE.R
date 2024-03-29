#' enrichGONE
#'
#' @description
#'  Prepares data obtained from DEG calculation to be suitable for Gene Enrichment Analysis using the enrichGO function of the clusterprofiler package.
#'  The two files are post filtered and log fold changed DEG data file which is input in the Universe parameter and fully DEG processed file including filtration, log fold change, padjustment and FDR correction.
#'  The latter file is used in the 'genes' attribute of enrichGO. The enrichGO object is used to generate the barplot, dotplot and upset plot for enrichment visualisation.
#'
#' @param input_file_1 This is the path to the file containing the differentially expressed genes (DEGs) for the genes input in enrichGO.
#' The data in this file should be pre-processed to filter out low-expressed genes and adjust for log fold change.

#' @param input_file_2 This is the path to the file containing the differentially expressed genes (DEGs) for the universe input in enrichGO.
#' The data in this file should be pre-processed to filter out low-expressed genes, adjust for log fold change, p-value, and FDR.
#'
#' @param output_file This is the path to the .xlsx file that will be created, containing the results of enrichGO analysis..
#'
#' @param
#'
#' @examples
#' library(DEGpathwayAnalysis)
#'
#'
#'
#'
#'

enrichGONE<- function(input_file_1,input_file_2,output_file) {


  library(clusterProfiler)
  library(edgeR)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
  library(ggupset)

#Importing file for gene attribute
#input_file_1<- "processed_data/DEGs_from_E-MTAB-2523.xlsx"
df = read.xlsx(input_file_1)


#Importing file for universe attribute
#input_file_2<- "processed_data/universe.xlsx"
df_universe = read.xlsx(input_file_2)



#FOR genes parameter in enrichgo-----------------------------------------

gene_list <- df$logFC


#Naming the genes based on the SYMBOL ID  found in the dataframe
names(gene_list) <- df[,1]

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# convert gene symbols to ENTREZID using bitr function
gene_list_entrez<- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db",drop=TRUE)

# check the converted gene list
#gene_list_entrez


#for universe parameter in enrichgo---------------------------------------


universe_gene_list = df_universe$logFC

#Naming the genes based on the SYMBOL ID  found in the dataframe

names(universe_gene_list) = df_universe[,1]


# sort the list in decreasing order (required for clusterProfiler)

universe_gene_list = sort(universe_gene_list,decreasing = TRUE)

#Conversion of SYMBOL ID to ENTREZID using bitr function

universe_list_entrez<- bitr(names(universe_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop=T)


go_enrich7 <- enrichGO(gene = gene_list_entrez$ENTREZID,universe = universe_list_entrez$ENTREZID, OrgDb = "org.Hs.eg.db" , keyType = 'ENTREZID',ont = "BP",pvalueCutoff = 0.05, qvalueCutoff = 0.1)



# displays a bar plot of the data in the "go_enrich7" dataset
show(barplot(go_enrich7))

# displays a dot plot of the data in the "go_enrich7" dataset
show(dotplot(go_enrich7))

# displays an upset plot of the data in the "go_enrich7" dataset
show(upsetplot(go_enrich7))

# displays a heat map plot of the data in the "go_enrich7" dataset, with fold change values determined by the "gene_list" dataset
show(heatplot(go_enrich7, foldChange = gene_list))

# Use the write.xlsx function to export the enrichGO results to an Excel file
write.xlsx(go_enrich7, output_file)

}

#Test our function from the DEGpathwayAnalysis folder
enrichGONE("processed_data/DEGs_from_E-MTAB-2523.xlsx","processed_data/universe.xlsx","processed_data/enrichGONE_results.xlsx")






















