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
  make_option(c("-g","--group"),
              help="sample group information")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
otu.nor.fil0 = read.table(opt$inputFile,header=T,sep='\t')

otu.nor.fil = otu.nor.fil0[,-1]
rownames(otu.nor.fil) = otu.nor.fil0[,1]
#head(otu.nor.fil)
#sample cluster
distance <- dist(t(otu.nor.fil), method="euclidean")
clusters <- hclust(distance)
#pdf("2.pdf",height = 5, width = 5)

##pdf 
outpdf = paste(opt$outputFile,".pdf",sep="")
pdf(outpdf,width=6,height=6)

plot(clusters, cex = 0.75,
     main="Sample dedrogram tree", 
     xlab="distance", ylab="Height", sub="")
#dev.off()
if(opt$group != 'notProvided'){
  group = strsplit(opt$group,",")[[1]]
  suppressPackageStartupMessages(library(RColorBrewer))
  #handle colors
  n = length(unique(group))
  if(n == 1){
    colors = rep("grey",ncol(otu.nor.fil))
    cols = "grey"
  }else if(n == 2){
    colors = as.character(factor(group,levels=unique(group),
                                 labels=c("blue","red")))
    cols = c("blue","red")
  }else{
    colors = as.character(factor(group,levels=unique(group),
                                 labels=brewer.pal(length(unique(group)),"Set3")))
    cols = brewer.pal(length(unique(group)),"Set3")
  }
  
  legend("topright",fill=cols,border="grey",
         legend = unique(group),box.col = "white")
  
  ## add group color
  rect(0.75:(ncol(otu.nor.fil) - 0.25),
       clusters$height[ncol(otu.nor.fil) - 1] ,
       1.25:(ncol(otu.nor.fil) + 0.25),
       clusters$height[ncol(otu.nor.fil) - 1] + 0.5,
       col = colors[clusters$order])
}

dev.off()


##png 
pngfile = paste(opt$outputFile, ".png",sep="")
png(pngfile,width=5000,height = 5000,units = "px",bg="white",res=800)

plot(clusters, cex = 0.75,
     main="Sample dedrogram tree", 
     xlab="distance", ylab="Height", sub="")
#dev.off()
if(opt$group != 'notProvided'){
  group = strsplit(opt$group,",")[[1]]
  suppressPackageStartupMessages(library(RColorBrewer))
  #handle colors
  n = length(unique(group))
  if(n == 1){
    colors = rep("grey",ncol(otu.nor.fil))
    cols = "grey"
  }else if(n == 2){
    colors = as.character(factor(group,levels=unique(group),
                                 labels=c("blue","red")))
    cols = c("blue","red")
  }else{
    colors = as.character(factor(group,levels=unique(group),
                                 labels=brewer.pal(length(unique(group)),"Set3")))
    cols = brewer.pal(length(unique(group)),"Set3")
  }
  
  legend("topright",fill=cols,border="grey",
         legend = unique(group),box.col = "white")
  
  ## add group color
  rect(0.75:(ncol(otu.nor.fil) - 0.25),
       clusters$height[ncol(otu.nor.fil) - 1] ,
       1.25:(ncol(otu.nor.fil) + 0.25),
       clusters$height[ncol(otu.nor.fil) - 1] + 0.5,
       col = colors[clusters$order])
}

dev.off()

