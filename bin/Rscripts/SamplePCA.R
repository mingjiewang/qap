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
  make_option(c("-d","--dimension"),
              help="PCA plot dimensions"),
  make_option(c("-s","--figStyle"),
              help="style of 2D scatter plot")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
otu.nor.fil0 = read.table(opt$inputFile,header=T,sep='\t')

otu.nor.fil = otu.nor.fil0[,-1]
rownames(otu.nor.fil) = otu.nor.fil0[,1]

#group and group color
group = strsplit(opt$group,",")[[1]]
#group = c(rep("case",5),rep("control",5))
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

pdffile = paste(opt$outputFile,".pdf",sep="")
pngfile = paste(opt$outputFile,".png",sep="")

##3d
if(opt$dimension == '3'){
  suppressPackageStartupMessages(library(scatterplot3d))
  
  if(nrow(otu.nor.fil) > ncol(otu.nor.fil)){
    pca = princomp(as.matrix(otu.nor.fil))
    plot.data <- data.frame(pca$loadings[,1:3])
  }else{
    fit <- prcomp(t(otu.nor.fil), scale=TRUE)
    plot.data <- data.frame(fit$x[,1:3])
  }
  colnames(plot.data) = c("PC1","PC2","PC3")
  
  diffangle <- function(ang){
    scatterplot3d(plot.data,main='PCA',color=colors,type='p',
                  highlight.3d=F,angle=ang,grid=T,box=T,scale.y=1,
                  cex.symbols=1.2,pch=16,col.grid='lightblue')
    legend("topright",
           unique(group),
           fill=cols,
           box.col="white")
  }
  
  #pdf file
  pdf(pdffile,onefile=TRUE,width=8,height=8)
  sapply(seq(-360,360,5),diffangle)
  dev.off()
  
  #png file
  png(pngfile, width=5000, height = 5000,
      units = "px", res=800, bg="white")
  diffangle(60)
  dev.off()
  
}

##2d
fit <- prcomp(t(otu.nor.fil), scale=TRUE)
plot.data <- data.frame(fit$x[,1:2])
suppressPackageStartupMessages(library(ggplot2))
if(opt$dimension == '2'){
  if(opt$figStyle == '1'){
    ## no circle, no label
    suppressPackageStartupMessages(library(ggplot2))
    p <- ggplot(plot.data,aes(PC1,PC2,shape=factor(group))) + 
      geom_point(aes(colour=factor(group)),size=2, stat="identity") +
      theme(panel.background = element_blank(),
            panel.border = element_rect(colour="black",fill=NA,size=0.5), 
            panel.grid.minor = element_line(colour = "#cccccc", size=0.1), 
            panel.grid.major = element_line(colour = "#cccccc", size=0.1)) +
      labs(title="PCA plot") +
      xlab("PC1") + ylab("PC2")
    
    ##save all plots
    p
    ggsave(pdffile,device="pdf",width=6,height=5,units = "in")
    p
    ggsave(pngfile,device="png",width=6,height=5,units = "in")
    
  }else if(opt$figStyle == '2'){
    ## with circle, no label
    suppressPackageStartupMessages(library(ggbiplot))
    p <- ggbiplot(fit, choices=1:2, 
                  var.axes=FALSE, 
                  obs.scale=1, var.scale=1, 
                  groups=group, 
                  ellipse=TRUE, circle=TRUE)
    
    #save all figures
    p
    ggsave(pdffile,device="pdf",width=5,height=6,units = "in")
    p
    ggsave(pngfile,device="png",width=5,height=6,units = "in")
    
  }else{
    ## no circle, with label
    suppressPackageStartupMessages(library(ggpubr))
    
    pcavalues <- princomp(otu.nor.fil)
    
    PCA2D_data <- data.frame(pcavalues$loadings[,1:2],
                             sample = colnames(otu.nor.fil),
                             group = factor(group))
    
    #pdf file
    p <- ggscatter(PCA2D_data, x = "Comp.1", y = "Comp.2",
              xlab = "PC1", ylab = "PC2",
              color = "group", shape = "group", 
              label = "sample",
              palette = cols,
              font.label = 8, main = "PCA", repel = FALSE,
              rug = TRUE, legend = "right")
    
    #save figures
    p
    ggsave(pdffile,device="pdf",width=6,height=5,units = "in")
    p
    ggsave(pngfile,device="png",width=6,height=5,units = "in")

  }
}
