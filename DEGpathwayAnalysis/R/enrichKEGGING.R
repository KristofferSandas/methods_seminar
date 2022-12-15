#' enrichKEGG
#'
#' @description
#'  Prepares data obtained from DEG calculation to be suitable for Gene Enrichment Analysis using the enrichKEGG function of the clusterprofiler package. 
#'  The two files are post filtered and log fold changed DEG data file which is input in the Universe parameter and fully DEG processed file including filtration, log fold change, padjustment and FDR correction. 
#'  The latter file is used in the 'genes' attribute of enrichKEGG.
#'
#' @param input_file_1 This is the path to the file containing the differentially expressed genes (DEGs) for the genes input in enrichKEGG. 
#' The data in this file should be pre-processed to filter out low-expressed genes and adjust for log fold change.

#' @param input_file_2 This is the path to the file containing the differentially expressed genes (DEGs) for the universe input in enrichKEGG. 
#' The data in this file should be pre-processed to filter out low-expressed genes, adjust for log fold change, p-value, and FDR.
#'
#' @param output_file This is the path to the .xlsx file that will be created, containing the results of enrichKEGG analysis..
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

enrichKEGGING<- function(input_file_1,input_file_2,output_file) {
  
  
  library(clusterProfiler)
  library(edgeR)
  library(openxlsx)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggnewscale)
  library(ggupset)
  
  #Importing file for gene attribute
  
  df = read.xlsx(input_file_1)
  
  
  #Importing file for universe attribute
  
  df_universe = read.xlsx(input_file_2)
  
  
  
  #FOR genes parameter in enrichKEGG-----------------------------------------
  
  gene_list <- df$logFC
  
  
  #Naming the genes based on the SYMBOL ID  found in the dataframe
  names(gene_list) <- df[,1]
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  # convert gene symbols to ENTREZID using bitr function
  gene_list_entrez<- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop=T)
  
  # check the converted gene list
  gene_list_entrez
  
  
  #for universe parameter in enrichKEGG---------------------------------------
  
  
  universe_gene_list = df_universe$logFC
  
  #Naming the genes based on the SYMBOL ID  found in the dataframe
  
  names(universe_gene_list) = df_universe[,1]
  
  
  # sort the list in decreasing order (required for clusterProfiler)
  
  universe_gene_list = sort(universe_gene_list,decreasing = TRUE)
  
  #Conversion of SYMBOL ID to ENTREZID using bitr function
  
  universe_list_entrez<- bitr(names(universe_gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db", drop=T)
  
  go_keggo <- enrichKEGG(gene = gene_list_entrez$ENTREZID,universe = universe_list_entrez$ENTREZID, organism  = "hsa", keyType = 'ncbi-geneid',pvalueCutoff = 0.05, qvalueCutoff = 0.1)
  

  
  # displays a bar plot of the data in the "go_keggo" dataset
  show(barplot(go_keggo))
  
  # displays a dot plot of the data in the "go_keggo" dataset
  show(dotplot(go_keggo))
  
  # displays an upset plot of the data in the "go_keggo" dataset
  show(upsetplot(go_keggo))
  
  # displays a heat map plot of the data in the "go_keggo" dataset, with fold change values determined by the "gene_list" dataset
  show(heatplot(go_keggo, foldChange = gene_list))
  
  
  # Use the write.xlsx function to export the enrichKEGG results to an Excel file
  write.xlsx(go_keggo, output_file)
  
}

#Test our function from the DEGpathwayAnalysis folder
enrichKEGGING("DEGs_from_E-MTAB-2523.xlsx","universe2.xlsx","kegging.xlsx")




