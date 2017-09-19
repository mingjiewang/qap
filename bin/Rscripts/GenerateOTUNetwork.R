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
              help="output result file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
otu.nor.fil0 = read.table(opt$inputFile,header=T,sep='\t')
#otu.nor.fil0 = read.table("NormalizedOTUCounts.txt",header=T,sep="\t")
otu.nor.fil = otu.nor.fil0[,-1]
rownames(otu.nor.fil) = otu.nor.fil0[,1]

#start to parse otu table
##----links
link.out = NULL
for (i in c(1:nrow(otu.nor.fil))){
  SourceNode.tmp = colnames(otu.nor.fil)[otu.nor.fil[i,] > 0]
  TargetNode.tmp = rep(rownames(otu.nor.fil)[i],length(SourceNode.tmp))
  Weight.tmp = unlist(otu.nor.fil[i,otu.nor.fil[i,] > 0])
  names(Weight.tmp) = NULL
  
  link.out$SourceNode = c(link.out$SourceNode,SourceNode.tmp)
  link.out$TargetNode = c(link.out$TargetNode,TargetNode.tmp)
  link.out$Weight = c(link.out$Weight,Weight.tmp)
}
link.out = as.data.frame(link.out)

write.table(link.out,paste(opt$outputFile,".Links.txt",sep=""),
            col.names = T,row.names = F,quote=F,sep="\t")

##--nodes
node.out = NULL
node.out$Node = c(colnames(otu.nor.fil),rownames(otu.nor.fil))
node.out$Attr.Category = c(rep("Sample",ncol(otu.nor.fil)),
                           rep("OTU",nrow(otu.nor.fil)))
node.out$Attr.Degree = c(colSums(otu.nor.fil > 0),
                         rowSums(otu.nor.fil > 0))
node.out = as.data.frame(node.out)

write.table(node.out,paste(opt$outputFile,".Nodes.txt",sep=""),
            col.names = T,row.names = F,quote=F,sep="\t")


