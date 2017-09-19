#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(gplots))
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-1","--inputFile1"),
              help="input nwk file"),
  make_option(c("-2","--inputFile2"),
              help="input consensus nwk file"),
  make_option(c("-o","--outputFile"),
              help="output result file"),
  make_option(c("-c","--colorScale"),
              help="the color of heatmap"),
  make_option(c("-f","--fontScale"),
              help="the degree of heatmap color"),
  make_option(c("-e","--edgeLabel"),
              help="add edge label or not"),
  make_option(c("-n","--nodeLabel"),
              help="add node label or not"),
  make_option(c("-a","--abundFile"),
              help="abundance file"),
  make_option(c("-s","--treeStyle"),
              help="style of tree plot")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
#nwk.file = file("../phylotest/res/tree.nwk","r")
nwk.file = file(opt$inputFile1,"r")
nwk.line = readLines(nwk.file, 1)
nwk.tree = read.tree(text = nwk.line)

#nwk.cs.file = file("../phylotest/res/tree_consensus.nwk","r")
nwk.cs.file = file(opt$inputFile2,"r")
nwk.cs.line = readLines(nwk.cs.file, 1)
nwk.cs.tree = read.tree(text = nwk.cs.line)

#dat.count = read.table("../phylotest/res/test.abund",header=T,sep="\t")
dat.count = read.table(opt$abundFile,header=T,sep="\t")

##settings
#--colorScale
user.tip.color = "black"
if(as.character(opt$colorScale) == "Y"){
  colors = colorpanel(length(nwk.tree$tip.label),"blue","red")
  count = dat.count$count[match(nwk.tree$tip.label,dat.count$seqID)]
  user.tip.color = colors[floor(rank(count))]
}

#--fontScale
user.font.size = 1
if(as.character(opt$fontScale) == "Y"){
  fonts = seq(0.4,1.6,length.out=length(nwk.tree$tip.label))
  count = dat.count$count[match(nwk.tree$tip.label,dat.count$seqID)]
  user.font.size = fonts[floor(rank(count))]
}

if(min(dat.count$count) == max(dat.count$count)){
  user.tip.color = "black"
  user.font.size = 1
}

##pdf tree
pdffile = paste(opt$outputFile,".pdf",sep="")
pdf(pdffile,width=7.5,height=5)
nf <- layout(matrix(c(1,1,1,1,1,1,1,1,2),3,3,byrow=TRUE), c(1.5,1.5), c(1.5,3))  
#layout.show(nf)
plot(nwk.tree, 
     type=as.character(opt$treeStyle),
     edge.width=1,
     label.offset = 0.2,
     font=2,
     tip.color = user.tip.color,
     cex=user.font.size
)

#--nodeLabel
if(as.character(opt$nodeLabel) == "Y"){
  node.labels = sprintf("%0.2f",as.numeric(nwk.tree$node.label))
  node.labels[node.labels == 'NA'] = ''
  nodelabels(node.labels,frame="rect",cex=0.6, bg="black", col="white", font=2)
}

#--edgeLabel
if(as.character(opt$edgeLabel) == "Y"){
  edge.labels = sprintf("%0.2f",as.numeric(nwk.tree$edge.length))
  edgelabels(edge.labels,frame="none",cex=0.5, col="black")
}

#figure legend
if(as.character(opt$colorScale) == "Y" & as.character(opt$fontScale) == "Y"){
  par(mar=c(3,2,2,3))
  barplot(c(100:300),space=0,
          col=colorpanel(201,"blue","red"),
          horiz=F,axes=F,border = NA,
          main="Sequence count",
          cex.main = 0.8)
  mtext("Low",2,las=1,cex=0.4,font=2)
  mtext("High",4,las=1,cex=0.8,font=2)
  axis(1,at = c(1,100,201),
       labels=c(min(dat.count$count),
                ceiling(mean(range(dat.count$count))),
                max(dat.count$count)),
       cex.axis=0.5) 
}
dev.off()

##pdf cs tree
pdffile = paste(opt$outputFile,"_Consensus.pdf",sep="")
pdf(pdffile,width=7.5,height=5)
nf <- layout(matrix(c(1,1,1,1,1,1,1,1,2),3,3,byrow=TRUE), c(1.5,1.5), c(1.5,3))  
#layout.show(nf)
plot(nwk.cs.tree, 
     type=as.character(opt$treeStyle),
     edge.width=1,
     label.offset = 0.2,
     font=2,
     tip.color = user.tip.color,
     cex=user.font.size
)

#--nodeLabel
if(as.character(opt$nodeLabel) == "Y"){
  node.labels = sprintf("%0.2f",as.numeric(nwk.cs.tree$node.label))
  node.labels[node.labels == 'NA'] = ''
  nodelabels(node.labels,frame="rect",cex=0.6, bg="black", col="white", font=2)
}

#figure legend
if(as.character(opt$colorScale) == "Y" & as.character(opt$fontScale) == "Y"){
  par(mar=c(3,2,2,3))
  barplot(c(100:300),space=0,
          col=colorpanel(201,"blue","red"),
          horiz=F,axes=F,border = NA,
          main="Sequence count",
          cex.main = 0.8)
  mtext("Low",2,las=1,cex=0.4,font=2)
  mtext("High",4,las=1,cex=0.8,font=2)
  axis(1,at = c(1,100,201),
       labels=c(min(dat.count$count),
                ceiling(mean(range(dat.count$count))),
                max(dat.count$count)),
       cex.axis=0.5) 
}
dev.off()

##png tree
pngfile = paste(opt$outputFile, ".png",sep="")
png(pngfile,width=7500,height = 5000,units = "px",bg="white",res=800)
nf <- layout(matrix(c(1,1,1,1,1,1,1,1,2),3,3,byrow=TRUE), c(1,1.5), c(1,3))  
layout.show(nf)
plot(nwk.tree, 
     type=as.character(opt$treeStyle),
     edge.width=1,
     label.offset = 0.2,
     font=2,
     tip.color = user.tip.color,
     cex=user.font.size
)

#--nodeLabel
if(as.character(opt$nodeLabel) == "Y"){
  node.labels = sprintf("%0.2f",as.numeric(nwk.tree$node.label))
  node.labels[node.labels == 'NA'] = ''
  nodelabels(node.labels,frame="rect",cex=0.6, bg="black", col="white", font=2)
}

#--edgeLabel
if(as.character(opt$edgeLabel) == "Y"){
  edge.labels = sprintf("%0.2f",as.numeric(nwk.tree$edge.length))
  edgelabels(edge.labels,frame="none",cex=0.5, col="black")
}

#figure legend
if(as.character(opt$colorScale) == "Y" & as.character(opt$fontScale) == "Y"){
  par(mar=c(3,2,2,3))
  barplot(c(100:300),space=0,
          col=colorpanel(201,"blue","red"),
          horiz=F,axes=F,border = NA,
          main="Sequence count",
          cex.main = 0.8)
  mtext("Low",2,las=1,cex=0.4,font=2)
  mtext("High",4,las=1,cex=0.8,font=2)
  axis(1,at = c(1,100,201),
       labels=c(min(dat.count$count),
                ceiling(mean(range(dat.count$count))),
                max(dat.count$count)),
       cex.axis=0.5) 
}
dev.off()

##png cs tree
pngfile = paste(opt$outputFile, "_Consensus.png",sep="")
png(pngfile,width=7500,height = 5000,units = "px",bg="white",res=800)
nf <- layout(matrix(c(1,1,1,1,1,1,1,1,2),3,3,byrow=TRUE), c(1,1.5), c(1,3))  
#layout.show(nf)
plot(nwk.cs.tree, 
     type=as.character(opt$treeStyle),
     edge.width=1,
     label.offset = 0.2,
     font=2,
     tip.color = user.tip.color,
     cex=user.font.size
)

#--nodeLabel
if(as.character(opt$nodeLabel) == "Y"){
  node.labels = sprintf("%0.2f",as.numeric(nwk.cs.tree$node.label))
  node.labels[node.labels == 'NA'] = ''
  nodelabels(node.labels,frame="rect",cex=0.6, bg="black", col="white", font=2)
}

#figure legend
if(as.character(opt$colorScale) == "Y" & as.character(opt$fontScale) == "Y"){
  par(mar=c(3,2,2,3))
  barplot(c(100:300),space=0,
          col=colorpanel(201,"blue","red"),
          horiz=F,axes=F,border = NA,
          main="Sequence count",
          cex.main = 0.8)
  mtext("Low",2,las=1,cex=0.4,font=2)
  mtext("High",4,las=1,cex=0.8,font=2)
  axis(1,at = c(1,100,201),
       labels=c(min(dat.count$count),
                ceiling(mean(range(dat.count$count))),
                max(dat.count$count)),
       cex.axis=0.5) 
}
dev.off()


