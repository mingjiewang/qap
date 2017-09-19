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
              help="sample group information"),
  make_option(c("-c","--heatColor"),
              help="the color of heatmap"),
  make_option(c("-d","--colorDegree"),
              help="the degree of heatmap color"),
  make_option(c("-s","--script"),
              help="the path of GenerateHeatmapColor.R")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
otu.nor.fil0 = read.table(opt$inputFile,header=T,sep='\t')

otu.nor.fil = otu.nor.fil0[,-1]
rownames(otu.nor.fil) = otu.nor.fil0[,1]
#head(otu.nor.fil)
## group colors
group = strsplit(opt$group,",")[[1]]
suppressPackageStartupMessages(library(RColorBrewer))
colors = as.character(factor(group,levels=unique(group),
                             labels=c("blue","red")))
if(length(unique(group)) >= 3){
  colors = as.character(factor(group,levels=unique(group),
                               labels=brewer.pal(length(unique(group)),"Set3")))
}

cols = c("blue","red")
if(length(unique(group)) >= 3){
  cols = brewer.pal(length(unique(group)),"Set3")
}

##pdf 
outpdf = paste(opt$outputFile,".pdf",sep="")
pdf(outpdf,width=6,height=6)

suppressPackageStartupMessages(library(gplots))
suppressWarnings(source(as.character(opt$script)))
#pdf("4.pdf",height = 5, width = 5)
heatmap.2(as.matrix(otu.nor.fil),
          col = GenerateHeatmapColor(opt$heatColor,as.numeric(opt$colorDegree)), 
          hclust=function(x) hclust(x,method = "ward.D2"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          scale = "row",dendrogram = 'both',
          key = TRUE, symkey = FALSE, density.info = "none", 
          trace = "none", cexRow = 0.4, cexCol = 0.8,
          main = paste("Heatmap"),
          ColSideColors = colors
)
legend("bottomleft",fill=cols,border="grey",
       legend = unique(group),box.col = "grey")

dev.off()


##png 
pngfile = paste(opt$outputFile, ".png",sep="")
png(pngfile,width=5000,height = 5000,units = "px",bg="white",res=800)

heatmap.2(as.matrix(otu.nor.fil),
          col = GenerateHeatmapColor(opt$heatColor,as.numeric(opt$colorDegree)), 
          hclust=function(x) hclust(x,method = "ward.D2"),
          distfun=function(x) as.dist((1-cor(t(x)))/2),
          scale = "row",dendrogram = 'both',
          key = TRUE, symkey = FALSE, density.info = "none", 
          trace = "none", cexRow = 0.4, cexCol = 0.8,
          main = paste("Heatmap"),
          ColSideColors = colors
)
legend("bottomleft",fill=cols,border="grey",
       legend = unique(group),box.col = "white")

dev.off()