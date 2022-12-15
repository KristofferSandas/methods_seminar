library(clusterProfiler)
library(openxlsx)
library(enrichplot)
library(DOSE)

# read in the results from enrichGO and enrichKEGG
go_results <- read.xlsx("enrichGONE_results.xlsx")
kegg_results <- read.xlsx("kegg_results.xlsx", 1)



# Create cneplot using 'ID', 'GeneRatio', and 'BgRatio' columns
# Set showCategory to 5, layout to "kk", and colorEdge to TRUE

enrichGO <- readEnrichment(go_results)

cnetplot(go_results, showCategory = 5, categorySize = "geneNum", foldChange = 2, fixed = TRUE)

