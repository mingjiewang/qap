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
              help="output result file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in data
dat = read.table(opt$inputFile,header=T,sep="\t")

#manipulate color
lblue = "#74d7ff"
red2 = "#ff1751"
grey = '#565656'

#manipulate data
dat$Interval = paste(dat$Start,dat$End,sep='-')
dat.heatmap = as.data.frame(dat[,c('MFI')])
rownames(dat.heatmap) = dat$Interval
colnames(dat.heatmap) = 'MFI'
dat.heatmap.t = t(dat.heatmap)
x.len=ncol(dat.heatmap.t)
y.len=nrow(dat.heatmap.t)
colors = colorpanel(x.len,low=lblue,high=red2)

##pdf file
#figure position 
pdffile = paste(opt$outputFile,".pdf",sep="",collapse = '')
pdf(pdffile,width = 20,height = 10)
nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(20,20), c(7,3), TRUE)  
#layout.show(nf)

#line plot
par(mar=c(5,5,3,1))
plot(rowMeans(dat[,c('Start','End')]),dat$MFI,xlab="Position",ylab="MFI",
     type="l",col=grey,pch=16,xaxt='n', main="Mutation Frequency Index",lwd=1.8)
points(rowMeans(dat[,c('Start','End')]),dat$MFI,col=colors[rank(dat$MFI)],pch=16,cex=1.5)
axis(1, at = c(1,dat$Start[2:nrow(dat)],dat$End[nrow(dat)]),
     labels = c(1, dat$Start[2:nrow(dat)] - 1,dat$End[nrow(dat)]))


#heatmap
par(mar=c(5,7,0,3))
max = max(dat.heatmap$MFI)
min = min(dat.heatmap$MFI)
image(x=1:x.len, y=1:y.len, t(dat.heatmap.t),col = colors,
      breaks=seq(min,max,by = (max - min) / x.len),
      axes = FALSE, ann=F)

abline(v=seq(0.5,0.5+x.len),col='white',lwd=0.25)
mtext("Heatmap",side = 2,las=1)
axis(1,seq(1,x.len),labels = dat$Interval,cex=0.2,las=2)
dev.off()

##png file
pngfile = paste(opt$outputFile,".png",sep="",collapse = '')
png(file = pngfile, width = 15000, height = 7500, res = 800)
nf <- layout(matrix(c(1,2),2,1,byrow=TRUE), c(20,20), c(7,3), TRUE)  
#layout.show(nf)

#line plot
par(mar=c(5,5,3,1))
plot(rowMeans(dat[,c('Start','End')]),dat$MFI,xlab="Position",ylab="MFI",
     type="l",col=grey,pch=16,xaxt='n', main="Mutation Frequency Index",lwd=1.8)
points(rowMeans(dat[,c('Start','End')]),dat$MFI,col=colors[rank(dat$MFI)],pch=16,cex=1.5)
axis(1, at = c(1,dat$Start[2:nrow(dat)],dat$End[nrow(dat)]),
     labels = c(1, dat$Start[2:nrow(dat)] - 1,dat$End[nrow(dat)]))


#heatmap
par(mar=c(5,7,0,3))
max = max(dat.heatmap$MFI)
min = min(dat.heatmap$MFI)
image(x=1:x.len, y=1:y.len, t(dat.heatmap.t),col = colors,
      breaks=seq(min,max,by = (max - min) / x.len),
      axes = FALSE, ann=F)

abline(v=seq(0.5,0.5+x.len),col='white',lwd=0.25)
mtext("Heatmap",side = 2,las=1)
axis(1,seq(1,x.len),labels = dat$Interval,cex=0.2,las=2)
dev.off()
