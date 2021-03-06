#load packages
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(RColorBrewer))
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
  make_option(c("-o","--outputFile"),
              help="output result file"),
  make_option(c("-p","--photos"),
              help="whether draw graphs"),
  make_option(c("-n","--selectSeqNum"),
              help="number of selected sequences"),
  make_option(c("-r","--repNum"),
              help="number of replicates"),
  make_option(c("-m","--normalized"),
              help="calculate normalized efficiency")
)
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
opt <- parse_args(OptionParser(option_list=option_list))


##read in file
dat0 = read.table(as.character(opt$inputFile), sep="\t", header=F)
#dat0 = read.table("1.txt",sep="\t",header = F)
dat = dat0$V1

##calculate
calculateSnEff <- function(dat, selectNum, repNum){
  if(selectNum < 1){
    num = length(dat)
    count = table(dat)
    freq = count / num
    freq.log = log(freq)
    shannon = sum(freq * freq.log)/log(num) * -1
    return(shannon)
  }else{
    out = NULL
    for (i in 1:repNum){
      if(selectNum > length(dat)){
        selectNum = length(dat)
      }
      dat.new = dat[sample(1:length(dat), selectNum, replace=F)]
      #print(paste(dat.new))
      num = length(dat.new)
      count = table(dat.new)
      freq = count / num
      freq.log = log(freq)
      shannon = sum(freq * freq.log)/log(num) * -1
      out = c(out,shannon)
    }
    return(mean(out))
  }
  
}


calculateSn <- function(dat){
  num = length(dat)
  count = table(dat)
  freq = count / num
  freq.log = log(freq)
  shannon = sum(freq * freq.log) * -1
  return(shannon)
}

if(as.character(opt$normalized) == 'Y'){
  out = data.frame(sample=opt$sampleLabel, 
                   value=calculateSnEff(dat,
                                        as.numeric(as.character(opt$selectSeqNum)), 
                                        as.numeric(as.character(opt$repNum))))
}else{
  out = data.frame(sample=opt$sampleLabel, 
                   value=calculateSn(dat))
}

##output
write.table(out,opt$outputFile,col.names=F,row.names=F,quote=F,sep="\t",append = T)

count = table(dat)
if(opt$photos == '1'){
  gra = matrix(sort(as.numeric(count),decreasing = T),ncol=1,byrow = F)
  
  if (length(count) >= 12){
    color = brewer.pal(12, 'Set3')
  }else if (length(count) >= 3){
    color = brewer.pal(length(count),'Set3')
  }else{
    color = brewer.pal(3,'Set3')
  }
  
  #pdf file
  pdf(file = paste(opt$sampleLabel,".pdf",collapse = '',sep=''),width = 10, height = 4)
  barplot(gra,horiz = T,beside = F, axes = F, col = color, border = F, 
          main=paste("Population Structure for", opt$sampleLabel))
  dev.off()
  #png file
  png(filename = paste(opt$sampleLabel, ".png", collapse = '',sep=''), width = 5000, height = 2000, res = 500)
  barplot(gra,horiz = T,beside = F, axes = F, col = color, border = F, 
          main=paste("Population Structure for", opt$sampleLabel))
  dev.off()
}else{
  #nothing
}



