source("https://bioconductor.org/biocLite.R")
options(repos = c(CRAN="http://cran.r-project.org"))
biocLite("Biostrings", ask = FALSE)
install.packages(c("gplots", "RColorBrewer", "ggplot2", "seqinr", "optparse", "plyr", 
                 "xlsx", "stringr", "scatterplot3d", "ggpubr"))
