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
  make_option(c("-w","--barWidth"),
              help="the width of bars"),
  make_option(c("-l","--legend"),
              help="whether show figure legend")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
otu.nor.fil0 = read.table(opt$inputFile,header=T,sep='\t')
#otu.nor.fil0 = read.table("NormalizedOTUCounts.txt",header = T,sep="\t")
otu.nor.fil = otu.nor.fil0[,-1]
rownames(otu.nor.fil) = otu.nor.fil0[,1]

library(reshape2)
dat = suppressMessages(melt(otu.nor.fil))
dat$otu = rep(rownames(otu.nor.fil),ncol(otu.nor.fil))
dat = dat[dat$value > 0,]
library(ggplot2)
p <- ggplot(dat, aes(variable, value))+
  geom_bar(aes(fill=otu),position="fill",
           stat="identity",
           width = as.numeric(opt$barWidth))+
  theme(axis.text.x=element_text(angle=45,size=10,
                                 face="bold",vjust=0.9,
                                 hjust=1),
        panel.background = element_blank(),
        panel.border = element_rect(colour="black",fill=NA,size=0.75),
        panel.grid.minor = element_line(colour = "#cccccc", size=0.1), 
        panel.grid.major = element_line(colour = "#cccccc", size=0.1))+
  xlab("Samples")+
  ylab("Normalized relative abundance")+
  ggtitle("Multiple samples QS abundance plot")
  #geom_hline(aes(yintercept=c(0.5)),col="white",size=0.5)+
  #geom_hline(aes(yintercept=c(0.25)),col="white",size=0.5)
  
if(opt$legend == "F"){
  p <- p + guides(fill=FALSE)
}

#save all figures
pdffile = paste(opt$outputFile,".pdf",sep="")
pngfile = paste(opt$outputFile,".png",sep="")


ggsave(pdffile,device="pdf",width=6,height=6,units = "in")
ggsave(pngfile,device="png",width=6,height=6,units = "in")

