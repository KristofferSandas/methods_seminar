# this is how you run the function and import the package

library(DEGpathwayAnalysis)


calculateDEGs("input_data/E-MTAB-2523.counts.txt",
               "processed_data/DEGs_from_E-MTAB-2523.xlsx",
               c(1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,0,1))

# 0=control/healhty, 1= disease
length(c(1,0,1,1,0,1,0,1,0,1,1,1,0,1,1,1,0,1))


