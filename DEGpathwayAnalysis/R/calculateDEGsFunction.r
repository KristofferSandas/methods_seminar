#' This script imports & filters the data and exports it to a .csv file
#' Should we add a prompt for the user to specify the files used?
#' This documentation needs updating later

calculateDEGs<- function(input_file, output_file, sample_groups) {

  library(edgeR)
  library(openxlsx)

  # importing counts
  gene_counts<- read.table(input_file, header=T, as.is=T,
                           row.names = 1, sep="\t")

  # inspecting data
  #dim(gene_counts)
  #head(gene_counts, 5)

  # importing sample data REMOVE
  #sample_table<- read.table(sample_table, header = TRUE,
  #                          as.is = TRUE, sep="\t")

  # checking that the counts and sample table are in the same order
  #colnames(gene_counts)
  #sample_table$sample

  # creating list for status: normal or carcinoma REMOVE
  #status<- factor(sample_table$disease,
                 # levels = c("normal", "carcinoma"))

  # sample groups from user ADD MORE EXPLANATION
  sample_groupings<- factor(sample_groups)

  # create an edgeR list object for DEGs
  DEG <- DGEList(counts = gene_counts,
                 group = sample_groupings)

  # inspect
  #DEG$samples

  # filter low expression genes by edgeR's recommended method
  # using values of at least 10 counts in at least 3 samples
  filter <- filterByExpr(y = DEG)
  DEG_filtered <- DEG[filter, keep.lib.sizes=F]

  #numbr_of_discarded_genes<- dim(DEG)[1]-dim(DEG_filtered)[1]

  # calculating edgeR's scaling factors for normalization
  DEG_scaling <- calcNormFactors(object = DEG_filtered)

  # inspecting scaling factors out of curiosity.
  # the product of all scaling factors should be 1
  #scaling_factors<- DEG_scaling$samples$norm.factors
  #prod(scaling_factors)
  # and it is.

  # fit the negative biomial model & estimate dispersion
  # straight out of the edgeR manual
  # no idea how this works under the hood
  DEG_modelled <- estimateDisp(y = DEG_scaling)

  # Pairwise exact test used for p-values
  DEG_stats <- exactTest(object = DEG_modelled)

  # calculate FDR (adjusted p-values)
  # comparison: positive logFC means higher expression
  # in carcinoma compared to normal
  DEG_tops = topTags(object = DEG_stats, n = "Inf")
  #head(DEG_tops$table)

  # filtering by FDR and logFC
  # cutoff for FDR: 0.05
  # cutoff for logFC: 1.1 - 1.5 seems standard
  DEG_by_FDR<- DEG_tops[DEG_tops$table$FDR<0.05,]
  #tail(DEG_by_FDR)
  #head(DEG_by_FDR)
  #dim(DEG_by_FDR)

  DEG_by_FDR_logFC<- DEG_by_FDR[DEG_by_FDR$table$logFC<(-1.3)|DEG_by_FDR$table$logFC>1.3,]
  #tail(DEG_by_FDR_logFC)
  #head(DEG_by_FDR_logFC)
  #dim(DEG_by_FDR_logFC)

  #head(DEG_by_FDR_logFC$table)

  # export to excel file
  write.xlsx(DEG_by_FDR_logFC$table,
             output_file,
             rowNames = TRUE)
}

