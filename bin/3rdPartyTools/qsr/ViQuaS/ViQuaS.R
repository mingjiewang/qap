source("viquas_files/ViquasSource.R")

ip = commandArgs(trailingOnly = TRUE)

o = as.integer(ip[3])
r = as.integer(ip[4])
richnessEst = as.integer(ip[5])
diversityRegionLength = as.integer(ip[6])

o <- ifelse(!is.na(o),o,3)                                                                                                                                                                                    
r <- ifelse(!is.na(r),r,0.7)
richnessEst <- ifelse(!is.na(richnessEst),richnessEst,0)

refSeq = read.fasta(file=ip[1],as.string=T,forceDNAtolower=F,seqonly=T)
refLength  = nchar(refSeq)
diversityRegionLength <- ifelse(!is.na(diversityRegionLength),diversityRegionLength,refLength)

ViQuaS(referenceFile=ip[1],readBAMFile=ip[2],o=o,r=r,richnessEst=richnessEst,diversityRegionLength=diversityRegionLength)