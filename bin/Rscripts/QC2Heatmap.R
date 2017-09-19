#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(gplots))
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i","--input"),
              help="input data"),
  make_option(c("-o","--output"),
              help="output figure")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in data
dat = read.table(opt$input,header=T,sep="\t")
#View(dat)
dat1 = dat[,-1]
rownames(dat1) = dat[,1]
#View(dat1)
dat2 = matrix(factor(unlist(dat1),levels = c('PASS','WARN','FAIL'),labels = c(0,1,2)), byrow = F,
              ncol = ncol(dat1))
dat2 = matrix(as.numeric(unlist(dat2)),byrow = F, ncol=ncol(dat1))
rownames(dat2) = rownames(dat1)
colnames(dat2) = colnames(dat1)

#colors
green = "#4cb049"
yellow = "#f4ed35"
red = "#e41e25"
grey = '#565656'

x.len=ncol(dat2)
y.len=nrow(dat2)
colors = colorpanel(3,low=green, mid = yellow, high=red)

##png file
pngfile = paste(opt$output,".png",sep="",collapse = '')
png(file = pngfile, width = 10000, height = 3000, res = 800)

par(mar=c(15,15,0,0))
image(x=1:x.len, y=1:y.len, t(dat2),col = colors,
      breaks=seq(-0.5,2.5,by = 1),
      axes = FALSE, ann=F)
abline(v=seq(0.5,0.5+x.len),col='white',lwd=1)
abline(h=seq(1.5,1.5+y.len),col="white",lwd=1)
axis(1,seq(1,x.len),labels = colnames(dat1),cex=0.2,las=2,tick=F)
axis(2,seq(1,y.len),labels = rownames(dat1),cex=0.1,las=1,tick=F)

dev.off()


##pdf file
pdffile = paste(opt$output,".pdf",sep="",collapse = '')
pdf(pdffile,width = 18,height = 4)

par(mar=c(15,15,0,0))
image(x=1:x.len, y=1:y.len, t(dat2),col = colors,
      breaks=seq(-0.5,2.5,by = 1),
      axes = FALSE, ann=F)
abline(v=seq(0.5,0.5+x.len),col='white',lwd=1)
abline(h=seq(1.5,1.5+y.len),col="white",lwd=1)
axis(1,seq(1,x.len),labels = colnames(dat1),cex=0.2,las=2,tick=F)
axis(2,seq(1,y.len),labels = rownames(dat1),cex=0.1,las=1,tick=F)

dev.off()

