#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gplots))
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
  make_option(c("-c","--cutoff"),
              help="cutoff value of strain ratio")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in data
data = read.table(opt$inputFile,header=T,sep="\t")
#data = read.table("../doman/tmp/P01_180_dominantStrain.fasta.rInput",header=T,sep="\t")

cutoff = as.numeric(opt$cutoff) * 100
data$Ratio = (data$Count / sum(data$Count)) * 100

##plot pdf 
pdffile = paste(opt$outputFile,".pdf",sep="",collapse = '')
pdf(pdffile,width = 12,height = 6)

layout(matrix(c(1,2),nrow=1,byrow = T),widths = c(5,5),heights = c(5,5))
#layout.show(lay)
par(mar=c(0,3,2,3))
##pie 1
data.minor = data[data$Ratio < cutoff,]
data.major = data[data$Ratio >= cutoff,]
data.minor.major.summary = data.frame(Strain=c("MultipleStrainMinor","MultipleStrainMajor"),
                        Count=c(sum(data.minor$Count),sum(data.major$Count)),
                        Ratio=c(sum(data.minor$Ratio),sum(data.major$Ratio)))
data.minor.major.summary$Ratio = sprintf("%0.2f",data.minor.major.summary$Ratio)
data.minor.major.summary$Label = paste(c("Minor strains\n","Major strains\n"),data.minor.major.summary$Count,"(",data.minor.major.summary$Ratio,"%)",sep="")
pie(data.minor.major.summary$Count, 
    data.minor.major.summary$Label,
    radius=0.8, 
    clockwise=T, 
    init.angle=90, 
    density=NULL, 
    col=c("grey","red"), 
    lty=1,
    border = "white",
    cex=0.8)
title(paste("Ratio of Major/Minor virus strains\n(Cut off value ",sprintf("%0.2f",cutoff),"%)"),
      cex.main = 0.8, font.main= 2, col.main= "black")

#pie 2
data.fil = data.major
data.fil$Ratio = sprintf("%0.2f",data.fil$Ratio)
data.fil$Label = paste(data.fil$Count,"(",data.fil$Ratio,"%)",sep="")
data.fil$Label2 = data.fil$Label
data.fil$Label2[data.fil$Ratio < 10] = NA #strian with ratio < 10% have no label
pie(data.fil$Count, 
    data.fil$Label2,
    radius=0.8, 
    clockwise=T, 
    init.angle=90, 
    density=NULL, 
    col=rainbow(nrow(data.fil)), 
    lty=1,
    border = "white",
    cex=0.8)
title(paste("Ratio of dominant virus strains\n(Only major strains shown)"),
      cex.main = 0.8, font.main= 2, col.main= "black")

dev.off()


#plot png
pngfile = paste(opt$outputFile,".png",sep="",collapse = '')
png(file = pngfile, width = 6000, height = 3000, res = 800)

par(mar=c(0,3,2,3))
layout(matrix(c(1,2),nrow=1,byrow = T),widths = c(5,5),heights = c(5,5))
#layout.show(lay)

##pie 1
data.minor = data[data$Ratio < cutoff,]
data.major = data[data$Ratio >= cutoff,]
data.minor.major.summary = data.frame(Strain=c("MultipleStrainMinor","MultipleStrainMajor"),
                                      Count=c(sum(data.minor$Count),sum(data.major$Count)),
                                      Ratio=c(sum(data.minor$Ratio),sum(data.major$Ratio)))
data.minor.major.summary$Ratio = sprintf("%0.2f",data.minor.major.summary$Ratio)
data.minor.major.summary$Label = paste(c("Minor strains\n","Major strains\n"),data.minor.major.summary$Count,"(",data.minor.major.summary$Ratio,"%)",sep="")
pie(data.minor.major.summary$Count, 
    data.minor.major.summary$Label,
    radius=0.8, 
    clockwise=T, 
    init.angle=90, 
    density=NULL, 
    col=c("grey","red"), 
    lty=1,
    border = "white",
    cex=0.8)
title(paste("Ratio of Major/Minor virus strains\n(Cut off value ",sprintf("%0.2f",cutoff),"%)"),
      cex.main = 0.8, font.main= 2, col.main= "black")

#pie 2
data.fil = data.major
data.fil$Ratio = sprintf("%0.2f",data.fil$Ratio)
data.fil$Label = paste(data.fil$Count,"(",data.fil$Ratio,"%)",sep="")
data.fil$Label2 = data.fil$Label
data.fil$Label2[data.fil$Ratio < 10] = NA #strian with ratio < 10% have no label
pie(data.fil$Count, 
    data.fil$Label2,
    radius=0.8, 
    clockwise=T, 
    init.angle=90, 
    density=NULL, 
    col=rainbow(nrow(data.fil)), 
    lty=1,
    border = "white",
    cex=0.8)
title(paste("Ratio of dominant virus strains\n(Only major strains shown)"),
      cex.main = 0.8, font.main= 2, col.main= "black")

dev.off()



