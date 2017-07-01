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
              help="output result file"),
  make_option(c("-l","--fileLabel"),
              help="label of input file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
dat = read.table(opt$inputFile,header=F,sep='\t')

#change it to char matrix
mat = matrix(unlist(strsplit(dat$V1,'')),nrow=nrow(dat), byrow = T)
mat.t = t(mat)

#get cs
cs = paste(c(">",opt$fileLabel,"\n"),collapse = '',sep='')
for (i in 1:nrow(mat.t)){
  tmp = names(sort(table(mat.t[i,]),decreasing = T)[1])
  cs = c(cs,tmp)
}

#output
cs.string = paste(cs,sep='',collapse = '')
write.table(cs.string,as.character(opt$outputFile),col.names = F,row.names = F,quote=F,sep='')
