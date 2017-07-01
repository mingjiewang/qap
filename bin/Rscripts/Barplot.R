#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gplots))
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i","--inputFile"),
              help="input Raw Data"),
  make_option(c("-l","--sampleLabel"),
              help = "sample Label"),
  make_option(c("-o","--outputDir"),
              help="output result folder")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#data 
dat = read.table(opt$inputFile, header=T,sep="\t")

#config colors
color.raw = colorpanel(10,low="#2bbcef",high="#ef3f6d")
dat$Color = cut(dat$Shannon, 
                breaks = seq(0,max(dat$Shannon),length.out = 11),
                labels = color.raw)

#output file
#pdf
pdffile = paste(opt$outputDir,"/pdf/",opt$sampleLabel,".pdf",sep="",collapse = '')
pdf(file = pdffile,width = 30, height = 5)
ggplot(dat,aes(x=Start.Position,y=Shannon))+
  geom_bar(stat="identity",position="identity", fill=dat$Color)+
  scale_x_continuous(breaks=seq(0, max(dat$End.Position), 10))  
dev.off()

#tif
tiffile = paste(opt$outputDir,"/tif/",opt$sampleLabel,".tif",sep="",collapse = '')
tiff(file = tiffile, width = 6000, height = 1000, res = 800, compression = 'lzw')
ggplot(dat,aes(x=Start.Position,y=Shannon))+
  geom_bar(stat="identity",position="identity", fill=dat$Color)+
  scale_x_continuous(breaks=seq(0, max(dat$End.Position), 10))  
dev.off()

