#load packages
suppressPackageStartupMessages(library(optparse))
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
dat = read.table(opt$inputFile,header=F,sep="\t")
#dat = read.table('0.depth',header=F,sep="\t")
colnames(dat) = c('ref','pos','depth')
#head(dat)

##filter1
##the start position and end position must have a extremly high depth than its adjancency postion
##so get hint from the adjancency difference
dat$depth1 = dat$depth
dat$depth1[dat$depth1 == 0] = 1
dat$depth2 = dat$depth1[c(2:nrow(dat),nrow(dat))]
ratio1 = (dat$depth2 - dat$depth)/dat$depth1
#plot(dat$pos,ratio1)
ratio2 = (dat$depth - dat$depth2)/dat$depth2
#plot(dat$pos,ratio2)

pos11 = dat$pos[ratio1 > 10]
pos11 = pos11 + 1
pos12 = dat$pos[ratio2 > 10]

pos1 = c(pos11,pos12)

##filter2
##As the depth should be uniform along the whole genome.
##thus the gap region which have no mapped reads must have a very large difference from the max depth
##get hint from the ratio
posfil = dat$pos[max(dat$depth)/dat$depth > 10]
posfil2 = posfil[c(2:length(posfil),length(posfil))]
pos21 = posfil[which((posfil - posfil2) < -1) + 1]
pos21 = pos21 - 1
pos22 = posfil[(posfil - posfil2) < -1]
pos22 = pos22 + 1

pos2 = c(pos22,pos21)

##filter3
##use k-means to cluster postion regions into several clusters
##the expect value is 3 clusters, with 2 gap clusters at genome start and end.
#n = length(pos2) / 2 + 1
#cluster = kmeans(posfil,n,iter.max=1000,algorithm = 'Hartigan-Wong')$cluster

merge = c(pos1,pos2) 
out = merge[duplicated(merge)]

write.table(out,opt$outputFile,col.names = F,row.names = F,quote = F,sep="\t")