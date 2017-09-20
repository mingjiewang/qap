#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
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
  make_option(c("-r","--sampleRatio"),
              help="Discard strains shown in less than this sample ratio"),
  make_option(c("-c","--cutoff"),
              help="sequence count cutoff")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
dat = read.table(opt$inputFile,header=F,sep='\t')
#dat = read.table("MergeSeq.RInput",header=F,sep="\t")
tbl = as.data.frame(table(dat$V1,dat$V2))

otu = NULL
otu.nor = NULL
samples = unique(dat$V2)
for (i in samples){
  tmp = tbl$Freq[tbl$Var2 == i]
  tmp.nor = tbl$Freq[tbl$Var2 == i] / sum(tmp)
  otu = cbind(otu,tmp)
  otu.nor = cbind(otu.nor,tmp.nor)
}
rownames(otu) = paste("OTU",c(1:nrow(otu)),sep="")
colnames(otu) = samples
rownames(otu.nor) = paste("OTU",c(1:nrow(otu)),sep="")
colnames(otu.nor) = samples


otu.seq = NULL
otu.seq$OTU = paste("OTU",c(1:nrow(otu)),sep="")
otu.seq$Sequence = unique(tbl$Var1)
otu.seq = as.data.frame(otu.seq)

write.table(otu.seq,str_replace(opt$outputFile,"OTUTable.txt","OTUSeq.txt"),quote=F,row.names = F,col.names = T,sep="\t")

## filter zero
count.not.zero <- function(x) ( sum(x >= as.numeric(opt$cutoff)) )
otu.notZeroCount <- apply(otu,1,count.not.zero)  # 1 according row ;2 according col
index <- otu.notZeroCount >= ceiling(ncol(otu) * as.numeric(opt$sampleRatio))

#index <- otu.notZeroCount >= ceiling(ncol(otu) * 0.3)
otu.fil = otu[index,]
otu.fil = as.data.frame(otu.fil)
otu.fil$OTUName = rownames(otu.fil)
otu.fil = otu.fil[,c(ncol(otu.fil),c(1:(ncol(otu.fil) - 1)))]

write.table(otu.fil,opt$outputFile,quote = F,sep="\t",col.names=T,row.names = F)

## normalized otu table
otu.nor.fil = otu.nor[index,]

## data transformation
freq.min = min(otu.nor.fil[otu.nor.fil != 0])
coef = NULL
for (i in c(1:20)){
  if(freq.min * (10 ^ i) > 1){
    coef = 10 ^ i
    break
  }
}
otu.nor.fil = otu.nor.fil * coef
otu.nor.fil[otu.nor.fil == 0] = 1
otu.nor.fil = log(otu.nor.fil,2)
otu.nor.fil = as.data.frame(otu.nor.fil)
otu.nor.fil$OTUName = rownames(otu.nor.fil)
otu.nor.fil = otu.nor.fil[,c(ncol(otu.nor.fil),c(1:(ncol(otu.nor.fil) - 1)))]

write.table(otu.nor.fil, str_replace(opt$outputFile,"OTUTable.txt","NormalizedOTUCounts.txt"),quote=F,row.names=F,col.names=T,sep="\t")
