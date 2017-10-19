#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RColorBrewer))
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i","--inputFile"),
              help="input Raw Data"),
  make_option(c("-o","--outputFile"),
              help="output result file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
dat = read.table(opt$inputFile,header=F,sep='\t')

#change it to char matrix
mat = matrix(unlist(strsplit(dat$V1,'')),nrow=nrow(dat), byrow = T)
mat.t = t(mat)

#shannon function
##calculate
shannon <- function(data){
  num = length(data)
  count = table(data)
  freq = count / num
  freq.log = log(freq)
  shannon = sum(freq * freq.log)/log(num) * -1
  return(shannon)
}

##output
output = NULL
for (i in c(1:nrow(mat.t))){
  sha.value = shannon(mat.t[i,])
  sha = sprintf("%0.10f",sha.value)
  start = i
  end = i + 1
  
  output$Start.Position = c(output$Start.Position,start)
  output$End.Position = c(output$End.Position,end)
  output$Shannon = c(output$Shannon,sha)
}

output = as.data.frame(output)

write.table(output,opt$outputFile,col.names=T,row.names=F,quote=F,sep="\t")

