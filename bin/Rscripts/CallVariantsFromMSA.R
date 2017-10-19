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
  make_option(c("-r","--refSeqFile"),
              help="reference sequence file"),
  make_option(c("-s","--startPos"),
              help="the start pos of mapped seq"),
  make_option(c("-c","--cutoff"),
              help="the cut off value of mutation freq")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#get ref name
refname = str_replace(opt$refSeqFile,'.*\\/','')
refname = str_replace(refname,'_[0-9]+\\.txt$','')

#read in file
dat = read.table(opt$inputFile,header=F,sep='\t')
ref = read.table(opt$refSeqFile,header=F,sep="\t")

#change dat to char matrix
mat = matrix(unlist(strsplit(dat$V1,'')),nrow=nrow(dat), byrow = T)
seq.mat = t(mat)

#change ref to char matrix
mat2 = matrix(unlist(strsplit(ref$V1,'')),nrow=nrow(ref), byrow = T)
ref.mat = t(mat2)

#start calculate
out = NULL
for (i in c(1:nrow(seq.mat))){
  tbl = data.frame(table(seq.mat[i,seq.mat[i,] != ref.mat[i,]]))
  if (nrow(tbl) == 0){
    next
  }
  tbl$Pos = i
  tbl$Ref = ref.mat[i,]
  out = rbind(out,tbl)
}

out = as.data.frame(out)
out$Reference = refname
out$Position = out$Pos + as.numeric(opt$startPos) - 1
out$Ref.Base = out$Ref
out$Alt.Base = out$Var1
out$Depth = ncol(seq.mat)
out$Frequency = out$Freq / out$Depth
out$Count = out$Freq

out2 = out[,c(5:11)]

#filter bases on cutoff
out2.fil = out2[out2$Frequency >= as.numeric(opt$cutoff), ]

#output
write.xlsx2(out2.fil, opt$outputFile, sheetName = "Mutation Profile", col.names = T, row.names = F, quote=F, sep="\t")
outputfile = str_replace(as.character(opt$outputFile),'\\.xlsx$','.txt');
write.table(out2.fil, outputfile, col.names = T, row.names = F,quote = F, sep="\t")
