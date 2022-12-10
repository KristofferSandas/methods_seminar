#' calculateDEGs
#'
#' @description
#' Filters and normalizes RNA-seq expression data, calculates statistically significant DEGs and exports them to an excel file. The DEGs are filtered by FDR < 0.05 and by absolute logFC value > 1.3. The low expression filtering parameters are: genes with at least 10 counts in at least 3 samples are kept.
#'
#' @param input_file This is the path to the expression counts file.
#'
#' @param output_file This is the path to the .xlsx file that will be created, containing a list of DEGs.
#'
#' @param sample_groups This is a vector of ones and zeros, for examlpe c(1,0,0,0,1). This vector specifies which samlpes(columns) are control/healthy(0) and which are disease(1).
#'
#' @examples
#' library(DEGpathwayAnalysis)
#'
#' calculateDEGs("input_data/E-MTAB-2523.counts.txt",
#'               "processed_data/DEGs_from_E-MTAB-2523.xlsx",
#'               c(1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,0,1))
#'
#'
#'
#'


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

