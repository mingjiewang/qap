suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(futile.logger))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(xlsx))
options(stringsAsFactors=FALSE)
option_list <- list(
  make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
              help="Print extra output [default]"),
  make_option(c("-q", "--quietly"), action="store_false",
              dest="verbose", help="Print little output"),
  make_option(c("-i","--inputFiles"),
              help="input Raw Data"),
  make_option(c("-o","--outputFile"),
              help="output result file"),
  make_option(c("-x","--outputXLSX"),
              help="output result xlsx file")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))

#read in file
files = unlist(strsplit(opt$inputFiles,","))

#merge data
dat = NULL
for (f in files){
  dat.tmp = read.table(f,header = F,sep="\t")
  mutation = paste(dat.tmp$V1,dat.tmp$V2,dat.tmp$V3,dat.tmp$V4,sep="^")
  
  dat$Mutation = c(dat$Mutation,mutation)
  dat$Freq = c(dat$Freq,dat.tmp$V5)
  dat$Dep = c(dat$Dep,dat.tmp$V6)
}
dat = as.data.frame(dat)

#split apply combine
out = ddply(dat,.(Mutation),summarise,Frequency = mean(Freq),Depth = ceiling(mean(Dep)))

#output mutation freq
write.table(out,opt$outputFile,quote=F,sep="\t",col.names = T,row.names = F)

#output mutation table
mut.table = NULL
mut.table$Reference = matrix(unlist(strsplit(out$Mutation,'^',fixed = T)),nrow=nrow(out),byrow = T)[,1]
mut.table$Position = as.numeric(matrix(unlist(strsplit(out$Mutation,'^',fixed = T)),nrow=nrow(out),byrow = T)[,2])
mut.table$Ref.Base = matrix(unlist(strsplit(out$Mutation,'^',fixed = T)),nrow=nrow(out),byrow = T)[,3]
mut.table$Alt.Base = matrix(unlist(strsplit(out$Mutation,'^',fixed = T)),nrow=nrow(out),byrow = T)[,4]
mut.table$Depth = out$Depth
mut.table$Frequency = out$Frequency
mut.table$Count = ceiling(out$Depth * out$Frequency)
mut.table = as.data.frame(mut.table)

write.xlsx2(mut.table, opt$outputXLSX, sheetName = "Mutation Profile", col.names = T, row.names = F, quote=F, sep="\t")
resfile = str_replace(as.character(opt$outputXLSX),'\\.xlsx$','.txt');
write.table(mut.table, resfile, col.names = T, row.names = F, quote = F, sep="\t")


