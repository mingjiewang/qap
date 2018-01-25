suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(stringr))
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i","--inputRawData"),
              help="input Raw Data"),
  make_option(c("-d","--inputIDFile"),
              help = "input ID File"),
  make_option(c("-o","--outputFile"),
              help="output result file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

##read in file
flog.info(paste("[R OUTPUT] Reading in read ID file:", opt$inputIDFile))
id = read.csv2(opt$inputIDFile, header=F, sep="\t")
flog.info(paste("[R OUTPUT] Reading in sequence file:", opt$inputRawData))
seq = read.csv2(opt$inputRawData, header=T, sep="\t")

flog.info(paste("[R OUTPUT] Totally", nrow(id), "IDs and", ceiling(nrow(seq) / 4), "sequences are inputted."))

seq$ID2 = str_replace(seq$ID,' .*','')

out = seq[match(id$V1, seq$ID2),]

if (ncol(out) == 5){
  outString = paste(out[,1], out[,2], out[,3], out[,4], sep="\n");
}else{
  outString = paste(str_replace(out[,1],'^@','>'), out[,2], sep="\n")
}

write.table(outString, opt$outputFile, col.names = F, row.names = F, quote=F, sep="\t")

