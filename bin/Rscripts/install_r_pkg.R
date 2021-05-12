source("http://bioconductor.org/biocLite.R")
options(repos = c(CRAN="http://cran.r-project.org"))
setRepositories()
setRepositories(ind = c(1:9))
biocLite("Biostrings", ask = FALSE)
#install.packages(c("gplots", "RColorBrewer", "ggplot2", "seqinr", "optparse", "plyr", 
#                "xlsx", "stringr", "scatterplot3d", "ggpubr", "futile.logger"), dependencies = TRUE)
print("===Installing gplots===")
install.packages("gplots", dependencies = TRUE) 

print("===Installing RColorBrewer===")
install.packages("RColorBrewer", dependencies = TRUE) 

print("===Installing ggplot2===")
install.packages("ggplot2", dependencies = TRUE) 

print("===Installing seqinr===")
install.packages("seqinr", dependencies = TRUE) 

print("===Installing optparse===")
install.packages("optparse", dependencies = TRUE) 

print("===Installing plyr===")
install.packages("plyr", dependencies = TRUE) 

print("===Installing xlsx===")
install.packages("xlsx", dependencies = TRUE) 

print("===Installing stringr===")
install.packages("stringr", dependencies = TRUE)

print("===Installing scatterplot3d===")
install.packages("scatterplot3d", dependencies = TRUE) 

print("===Installing ggpubr===")
install.packages("ggpubr", dependencies = TRUE)

print("===Installing futile.logger===")
install.packages("futile.logger", dependencies = TRUE) 
