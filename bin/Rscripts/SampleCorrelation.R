#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(ape))
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
  make_option(c("-c","--heatColor"),
              help="color of heatmaps"),
  make_option(c("-b","--corTable"),
              help="whether output corr table"),
  make_option(c("-l","--withLabel"),
              help="whether add figure label"),
  make_option(c("-a","--triangle"),
              help="whether draw a triangle"),
  make_option(c("-s","--sort"),
              help="sort samples or not"),
  make_option(c("-d","--dendrogram"),
              help="whether draw dendrogram or not"),
  make_option(c("-e","--line"),
              help="draw white lines or not")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
otu.nor.fil0 = read.table(opt$inputFile,header=T,sep='\t')
#otu.nor.fil0 = read.table("NormalizedOTUCounts.txt",header=T,sep="\t")
otu.nor.fil = otu.nor.fil0[,-1]
rownames(otu.nor.fil) = otu.nor.fil0[,1]


##corr
cor = data.frame(Sample=colnames(otu.nor.fil))
for(i in c(1:ncol(otu.nor.fil))){
  cor.row = NULL
  for(k in c(1:ncol(otu.nor.fil))){
    #cor.tmp = NA
    #if(i == k){
      #cor.tmp = NA
    #}else{
      cor.tmp = cor(otu.nor.fil[,i],otu.nor.fil[,k])
    #}
    cor.row = c(cor.row,cor.tmp)
  }
  cor.row = as.data.frame(cor.row)
  cor = cbind(cor,cor.row)
}

cor2 = cor[,-1]
rownames(cor2) = cor[,1]
colnames(cor2) = cor[,1]

if(as.character(opt$sort) == 'Y'){
  distance <- dist(t(otu.nor.fil), method="euclidean")
  clusters <- hclust(distance)
  sample.order = clusters$order
  
  cor2 = cor2[sample.order,sample.order] 
}  
samples = colnames(cor2)
#cor2 = cor2[,c(ncol(cor2):1)]
#cor2[is.na(cor2)] = NULL

cor2.half = cor2
for (i in c(1:nrow(cor2))){
  for (k in c(1:nrow(cor2))){
      if(k > i){
        cor2.half[i,k] = NA
      }else{
        next
      }
    
  }
}


xlen=ncol(cor2)
ylen=nrow(cor2)
cor2.nona = cor2[!is.na(cor2)]
breaks = seq(min(cor2.nona),max(cor2.nona),(max(cor2.nona) - min(cor2.nona)) / (ncol(cor2) * 20))
cor2.label = matrix(sprintf("%0.2f",unlist(t(cor2))),nrow=nrow(cor2),byrow = F)

## outpdf 
pdffile = paste(opt$outputFile,".pdf",sep="")
pdf(pdffile,height=6,width = 6)


if(as.character(opt$triangle) == 'Y'){
  if(as.character(opt$dendrogram) == 'Y'){
    nf <- layout(matrix(c(1,2,2,2,3,4,4,4,3,4,4,4,3,4,4,4),4,4,byrow=TRUE), 
                 c(1,1), c(1,1), TRUE)  
    #layout.show(nf)
    
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",axes = F,xlab="",ylab="")
    
    plot(clusters, cex = 0.5, main = NULL, xlab = NULL, 
         axes = F, hang = -1, sub="", ylab="")
    
    plot(as.phylo(clusters), cex=0.5)
    
    par(mar=c(0.6,0.6,0.6,0.6))
    image(1:xlen,1:ylen,t(t(cor2.half[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
  }else{
    par(mar=c(1,5,5,1))
    image(1:xlen,1:ylen,t(t(cor2.half[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
    axis(3,c(1:ncol(cor2)),samples,tick = F,las=2,cex.axis=0.5)
    axis(2,c(1:nrow(cor2)),samples[nrow(cor2):1],tick = F,las=1,cex.axis=0.5)
  }
}else{
  if(as.character(opt$dendrogram) == 'Y'){
    nf <- layout(matrix(c(1,2,2,2,3,4,4,4,3,4,4,4,3,4,4,4),4,4,byrow=TRUE), 
                 c(1,1), c(1,1), TRUE)  
    #layout.show(nf)
    
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",axes = F,xlab="",ylab="")
    
    plot(clusters, cex = 0.5, main = NULL, xlab = NULL, 
         axes = F, hang = -1, sub="", ylab="")
    
    plot(as.phylo(clusters), cex=0.5)
    
    par(mar=c(0.6,0.6,0.6,0.6))
    image(1:xlen,1:ylen,t(t(cor2[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
  }else{
    par(mar=c(1,5,5,1))
    image(1:xlen,1:ylen,t(t(cor2[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
    axis(3,c(1:ncol(cor2)),samples,tick = F,las=2,cex.axis=0.5)
    axis(2,c(1:nrow(cor2)),samples[nrow(cor2):1],tick = F,las=1,cex.axis=0.5)
  }
}

if(as.character(opt$line) == 'Y'){
  abline(v = c(0.5:(ncol(cor2) + 0.5)),
         h=c(0.5:(ncol(cor2) + 0.5)),
         col="white")
}

if(opt$withLabel == 'Y'){
  text(x = rep(1:xlen,each=ncol(cor2)), 
       y = rep(ylen:1,ncol(cor2)),
       labels=unlist(t(cor2.label)), col="white",
       cex=0.5)
}

dev.off()

##png file
pngfile = paste(opt$outputFile,".png",sep="")
png(pngfile,width=5000,height=5000,res=800,units = 'px',bg="white")



if(as.character(opt$triangle) == 'Y'){
  if(as.character(opt$dendrogram) == 'Y'){
    nf <- layout(matrix(c(1,2,2,2,3,4,4,4,3,4,4,4,3,4,4,4),4,4,byrow=TRUE), 
                 c(1,1), c(1,1), TRUE)  
    #layout.show(nf)
    
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",axes = F,xlab="",ylab="")
    
    plot(clusters, cex = 0.5, main = NULL, xlab = NULL, 
         axes = F, hang = -1, sub="", ylab="")
    
    plot(as.phylo(clusters), cex=0.5)
    
    par(mar=c(0.6,0.6,0.6,0.6))
    image(1:xlen,1:ylen,t(t(cor2.half[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
  }else{
    par(mar=c(1,5,5,1))
    image(1:xlen,1:ylen,t(t(cor2.half[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
    axis(3,c(1:ncol(cor2)),samples,tick = F,las=2,cex.axis=0.5)
    axis(2,c(1:nrow(cor2)),samples[nrow(cor2):1],tick = F,las=1,cex.axis=0.5)
  }
}else{
  if(as.character(opt$dendrogram) == 'Y'){
    nf <- layout(matrix(c(1,2,2,2,3,4,4,4,3,4,4,4,3,4,4,4),4,4,byrow=TRUE), 
                 c(1,1), c(1,1), TRUE)  
    #layout.show(nf)
    
    par(mar=c(0,0,0,0))
    plot(1,1,type="n",axes = F,xlab="",ylab="")
    
    plot(clusters, cex = 0.5, main = NULL, xlab = NULL, 
         axes = F, hang = -1, sub="", ylab="")
    
    plot(as.phylo(clusters), cex=0.5)
    
    par(mar=c(0.6,0.6,0.6,0.6))
    image(1:xlen,1:ylen,t(t(cor2[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
  }else{
    par(mar=c(1,5,5,1))
    image(1:xlen,1:ylen,t(t(cor2[c(ncol(cor2):1),c(ncol(cor2):1)])),
          col=colorpanel(length(breaks) - 1,as.character(opt$heatColor),"orange","red"),
          #col=terrain.colors(length(breaks) - 1),
          breaks = breaks,
          axes = FALSE, ann=F)
    
    axis(3,c(1:ncol(cor2)),samples,tick = F,las=2,cex.axis=0.5)
    axis(2,c(1:nrow(cor2)),samples[nrow(cor2):1],tick = F,las=1,cex.axis=0.5)
  }
}

if(as.character(opt$line) == 'Y'){
  abline(v = c(0.5:(ncol(cor2) + 0.5)),
         h=c(0.5:(ncol(cor2) + 0.5)),
         col="white")
}

if(opt$withLabel == 'Y'){
  text(x = rep(1:xlen,each=ncol(cor2)), 
       y = rep(ylen:1,ncol(cor2)),
       labels=unlist(t(cor2.label)), col="white",
       cex=0.5)
}


dev.off()



##output corr table
if(opt$corTable == 'Y'){
  write.table(cor2,paste(opt$outputFile,".txt",sep=""),col.names = T,
              row.names = T,sep="\t",quote=F)
}
