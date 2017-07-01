#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(xlsx))
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
  make_option(c("-l","--seqLen"),
              help="length of sequence")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

dat = read.table(opt$inputFile, header=T,sep="\t")

numberOfMutations = nrow(dat)

mfi = numberOfMutations / (as.numeric(opt$seqLen) * mean(dat$Depth))

write.table(mfi,opt$outputFile,quote=F,sep="\t",col.names=F,row.names=F)

