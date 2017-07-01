library(Biostrings)
library(seqinr)

ViQuaS<-function(referenceFile,readBAMFile,o=3,r=0.7,richnessEst=0,diversityRegionLength)
{ 
  baseDir = getwd()
  
  path = "viquas";
  dir.create("Samples")
  dir.create("Samples/viquas")
  
  dir.create("SIM")
  
  refSeq = read.fasta(file=referenceFile,as.string=T,forceDNAtolower=F,seqonly=T)
  
  refLength  = nchar(refSeq) 
  
  write.fasta(sequences=refSeq,names="reference",file.out="SIM/reference.fsa",open="w")
  
  ri = alignFilter(baseDir, path, readBAMFile);
  
  readLength = ri[1]
  nReads = ri[2]
  
  #03
  olcCluster(baseDir, path, o, r);
  
  #04
  mutationCaller("SIM/reference.fsa", paste(path,"/contigs",path,".fsa",sep=""), path);
  
  #05
  st <- as.matrix(read.delim(file=paste(path,"/startpositions",path,".txt",sep=""), header=FALSE));
  cl <- as.matrix(read.delim(file=paste(path,"/contigs",path,"_lengths.txt",sep=""), header=FALSE));
  ln <- as.matrix(read.delim(file=paste(path,"/opLinks",path,".txt",sep=""), header=FALSE,sep=","));
  sn <- matrix(data=unlist(strsplit(ln[,2],split="-")),ncol=3,byrow=T);
  sn <- cbind(ln[,1],sn);
  nm <- as.matrix(unique(ln[,1]))
  tr <- as.matrix(read.delim(file=paste(path,"/contigs",path,"_readcounts.txt",sep=""), header=FALSE));
  co <- as.matrix(read.delim(file=paste(path,"/contigs",path,"_coverages.txt",sep=""), header=FALSE));  
  ct <- as.matrix(unlist(read.fasta(file=paste(path,"/contigs",path,".fsa",sep=""),as.string=T,forceDNAtolower=F,seqonly=T)));
  
  ctln <- as.matrix(read.delim(file=paste(path,"/opLinksCt",path,".txt",sep=""), header=FALSE,sep=","));
  ctsn <- matrix(data=unlist(strsplit(ctln[,2],split="-")),ncol=3,byrow=T);
  ctsn[ctsn[,3]=="I",1] <- as.integer((as.integer(ctsn[ctsn[,3]=="I",1])+(nchar(ctsn[ctsn[,3]=="I",2])-1)/2))+1
  ctsn <- cbind(ctln[,1],ctsn);  
  
  for(n in 1:nrow(nm))
  {
    snpList = matrix(data=unlist(sn[which(sn[,1]==nm[n,1]),]),ncol=4,byrow=F);
    ctMutationPositions = matrix(data=unlist(ctsn[which(ctsn[,1]==nm[n,1]),]),ncol=4,byrow=F);
    c = as.integer(substr(nm[n,1],start=7,stop=nchar(nm[n,1])));
    correctClusterError(paste(path,"/",nm[n,1],".txt",sep=""),st[c],cl[c],snpList,co[c],tr[c],path,ct[c],readLength,ctMutationPositions);
  }
  
  #06
  estimateFreq(path, paste(path,"/opclsters",path,".txt",sep=""),paste(path,"/commonReadAlnStEnd",path,".txt",sep=""),st);
  
  #07
  createMutationGraph(path,paste(path,"/opclsters",path,".txt",sep=""),paste(path,"/islandsFile",path,".txt",sep=""),paste(path,"/opclstersFrq",path,".txt",sep=""),paste(path,"/opclstersSeq",path,".fsa",sep=""),refLength);
  
  #08
  constructStrainsDynamic(path,1);
  
  #10
  if(as.integer(richnessEst) == 1)
  {
    QNE(path,refLength,readLength,nReads,diversityRegionLength)
  }
  
  system(paste("rm -rf Samples"))
  system(paste("rm -rf SIM"))
  system(paste("rm -rf viquas"))
}

#######################################################################################################

alignFilter<-function(baseDir, path, readBAMFile)
{
  system(paste("samtools view -h -o aln.sam ",readBAMFile,sep=""));
  #system(paste("samtools sort ",readBAMFile," alns"));
  system("samtools sort aln.sam > alns.bam");
  system("samtools index alns.bam");
  
  system(paste("cut -f 1 aln.sam >names.txt",sep=""));
  system(paste("cut -f 10 aln.sam >seq.txt",sep=""));
  df <- try(names <- read.table(file="names.txt",header=F,sep="",skip=2),silent=TRUE)
  if (class(df)=='try-error'){  names = 0 }
  df <- try(seq <- read.table(file="seq.txt",header=F,sep="",skip=2),silent=TRUE)
  if (class(df)=='try-error'){  seq = "A" }  
  write.table(x=names,file="names.txt",append=F,quote=F,row.names=F,col.names=F)
  write.table(x=seq,file="seq.txt",append=F,quote=F,row.names=F,col.names=F)  
  system(paste("viquas_files/forwardreads.pl ",path," names.txt seq.txt 1>/dev/null",sep=""));
  
  system(paste("viquas_files/exactfilter.pl -f ",path,"-readsnew.fa 1>/dev/null",sep=""));
  
  system("viquas_files/smalt_i686 index -k 8 -s 2 mmm SIM/reference.fsa 2>/dev/null");
  system(paste("viquas_files/smalt_i686 map -f samsoft -o alnc.sam mmm commonreads.fsa",sep=""));
  
  system(paste("cut -f 4 alnc.sam >alnpos.txt",sep=""));
  
  df <- try(al <- read.table(file="alnpos.txt",header=F,sep="",skip=2),silent=TRUE)
  if (class(df)=='try-error'){  al = 0 } 
  df <- try(rl <- read.table(file="rdlen.txt",header=F,sep=""),silent=TRUE)
  if (class(df)=='try-error'){  rl = 0 } 
  
  alnpos = as.matrix(al);
  alnrl  = as.matrix(rl);
  df <- try(alnend <- alnpos + alnrl,silent=TRUE)
  if (class(df)=='try-error'){  alnpos=0;alnend=0; }
  
  readLength = round(mean(nchar(as.character(seq$V1))),0)
  nReads = length(as.character(seq$V1))
  
  write.table(x=cbind(alnpos,alnend),file="commonReadAlnStEnd.txt",append=F,sep=",",row.names=F,col.names=F)
  
  system("rm mmm.smi 2>/dev/null");
  system("rm mmm.sma 2>/dev/null");
  system("rm aln.sam 2>/dev/null");
  system("rm alns.bam.bai 2>/dev/null");
  system("rm names.txt 2>/dev/null");
  system("rm seq.txt 2>/dev/null");
  system("rm rdlen.txt 2>/dev/null");
  system("rm alnc.sam 2>/dev/null");
  system(paste("rm ",path,"-readsnew.fa",sep=""));
  
  system(paste("mv filtered.fsa Samples/",path," 2>/dev/null",sep=""));
  system(paste("mv alns.bam Samples/",path," 2>/dev/null",sep=""));
  system(paste("mv alnpos.txt Samples/",path,"/alnpos",path,".txt"," 2>/dev/null",sep=""));
  system(paste("mv commonreads.fsa Samples/",path,"/commonreads",path,".fsa"," 2>/dev/null",sep=""));
  system(paste("mv commonReadAlnStEnd.txt Samples/",path," 2>/dev/null",sep=""));
  
  return(c(readLength,nReads))
}

olcCluster<-function(baseDir, path, o=3, r=0.7)
{
  system(paste("cp -R Samples/",path,"/filtered.fsa ",baseDir,"/",sep="")); 
  system(paste("viquas_files/SSAKE.pl -f filtered.fsa -w 1 -c 1 -z 1 -o ",as.character(o)," -r ",as.character(r)," 1>/dev/null",sep=""))
  system('viquas_files/edge.pl -f contigs.readposition 1>/dev/null')
  
  dir.create(path);
  system(paste("cp -R Samples/",path,"/commonReadAlnStEnd.txt ",path,"/commonReadAlnStEnd",path,".txt"," 2>/dev/null",sep=""));
  
  system(paste('mv filtered.fsa ',path,'/filtered',path,'.fsa',sep=''));
  system('rm contigs.log')
  system('rm contigs.short')
  system('rm contigs.singlets')
  system(paste('mv','contigs.fsa',paste(path,'/contigs',path,'.fsa',sep=""),sep=" "));
  system(paste('mv','contigs.coverage.csv',path,sep=" "));
  system(paste('mv','contigs.readposition',path,sep=" "));
  system(paste('mv','contig*.txt',path,sep=" "));
}

mutationCaller<-function(referenceFile, contigFile, path)
{  
  reference = readDNAStringSet(referenceFile, "fasta");
  contigs = readDNAStringSet(contigFile, "fasta");
  
  contigLengths = width(contigs);
  contigLengthsFile = paste(substr(contigFile,1,nchar(contigFile)-4),"_lengths.txt",sep="");
  write.table(x=as.matrix(contigLengths),file=contigLengthsFile,row.names=F,col.names=F,append=FALSE);
  
  contigDetails = matrix(unlist(strsplit(names(contigs),split="[|]")),nrow=length(contigLengths),byrow=T);
  contigNames = as.matrix(contigDetails[,1]);
  contigRds = as.matrix(as.numeric(substr(contigDetails[,3],5,width(contigDetails[,4]))));
  contigCov = as.matrix(as.numeric(substr(contigDetails[,4],4,width(contigDetails[,4]))));
  
  contigNamesFile = paste(substr(contigFile,1,nchar(contigFile)-4),"_names.txt",sep="");
  write.table(x=contigNames,file=contigNamesFile,row.names=F,col.names=F,append=FALSE);
  
  contigRdsFile = paste(substr(contigFile,1,nchar(contigFile)-4),"_readcounts.txt",sep="");
  write.table(x=contigRds,file=contigRdsFile,row.names=F,col.names=F,append=FALSE);
  
  contigCovFile = paste(substr(contigFile,1,nchar(contigFile)-4),"_coverages.txt",sep="");
  write.table(x=contigCov,file=contigCovFile,row.names=F,col.names=F,append=FALSE);
  
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE);
  
  for(i in 1:length(contigs))
  {
    localAlign = pairwiseAlignment(contigs[i], reference[1],type="overlap" ,substitutionMatrix = mat, gapOpening = -5, gapExtension = -2);
    
    snps = mismatchSummary(localAlign);
    st = as.integer(snps$subject$SubjectPosition);
    wd = c(rep(1,length(snps$subject$SubjectPosition)));
    ty = snps$subject$Pattern;
    opsnps = data.frame(st,wd,ty);
    st = as.integer(snps$subject$SubjectPosition) - start(subject(localAlign)) + start(pattern(localAlign));
    opsnpsCon = data.frame(st,wd,ty);
    
    loc = NULL;
    locCon = NULL;
    locRef = NULL;
    
    indels = indel(localAlign);
    ins = indels@insertion@unlistData;     
    del = indels@deletion@unlistData;
    
    if(length(ins)+length(del) > 0)
    {
      insRows = NULL;
      delRows = NULL;
      
      if(length(ins) > 0)
      {
        insRows = cbind(ins@start,ins@width, "I")
      }
      
      if(length(del) > 0)
      {
        delRows = cbind(del@start,del@width, "D")
      }
      
      locs = rbind(insRows,delRows)
      st = as.integer(locs[,1]);wd = as.integer(locs[,2]);ty = locs[,3];sq="";
      loc = data.frame(st,wd,ty,as.character(sq),stringsAsFactors=F)
      loc = loc[order(loc$st),]
      locRef = data.frame(st,wd,ty,as.character(sq),stringsAsFactors=F)
      locRef = locRef[order(locRef$st),]
      locCon = data.frame(st,wd,ty,as.character(sq),stringsAsFactors=F)
      locCon = locCon[order(locCon$st),]
      insoffset = 0;
      deloffset = 0;
      
      for(j in 1:nrow(locs))
      {
        if(j>1)
        {
          if(loc[j-1,3] == "I")
          {
            insoffset = insoffset + as.integer(loc[j-1,2]);
          }         
          if(locRef[j-1,3] == "D")
          {
            deloffset = deloffset + as.integer(loc[j-1,2]);
          }
        }         
        if(locRef[j,3] == "D")
        {
          locRef[j,1] = as.integer(locRef[j,1]) + start(subject(localAlign)) - start(pattern(localAlign)) - insoffset + deloffset - 1;
          locCon[j,1] = as.integer(locCon[j,1]) - 1;  
          
          mu = substr(x=reference[1],start=locRef[j,1]+1,stop=locRef[j,1]+as.integer(locRef[j,2]))
          locRef[j,4] = as.character(mu)
          locCon[j,4] = as.character(mu)
        } 
        if(locRef[j,3] == "I")
        {
          locRef[j,1] = as.integer(locRef[j,1]) + start(subject(localAlign)) - start(pattern(localAlign)) - 1;
          locCon[j,1] = as.integer(locCon[j,1]) + insoffset - deloffset - 1;
          
          mu = substr(x=contigs[i],start=locCon[j,1]+1,stop=locCon[j,1]+as.integer(locCon[j,2]))
          locRef[j,4] = as.character(mu)
          locCon[j,4] = as.character(mu)
        }
      }
    }
    
    st = locCon$st;wd = locCon$as.character.sq.;ty=locCon$ty;
    contigPositions = rbind(opsnpsCon,cbind(st,wd,ty));
    contigPositions = contigPositions[order(as.integer(contigPositions$st)),];
    
    st = locRef$st;wd = locRef$as.character.sq.;ty=locRef$ty;
    mutations = rbind(opsnps,cbind(st,wd,ty));
    mutations = mutations[order(as.integer(mutations$st)),];
    
    if(length(mutations$st)>0)
    {
      write(x=paste(contigPositions$st,contigPositions$wd,contigPositions$ty,sep="-"),file=paste(path,"/opMutationsCt",path,".txt",sep=""),ncolumns=1,append=T,sep=",");
      write.table(x=cbind(rep.int(x=contigNames[i],times=length(contigPositions$st)),paste(contigPositions$st,contigPositions$wd,contigPositions$ty,sep="-")),file=paste(path,"/opLinksCt",path,".txt",sep=""),append=T,sep=",",col.names=F,row.names=F,quote=F);
      write(x=paste(mutations$st,mutations$wd,mutations$ty,sep="-"),file=paste(path,"/opMutations",path,".txt",sep=""),ncolumns=1,append=T,sep=",");
      write.table(x=cbind(rep.int(x=contigNames[i],times=length(mutations$st)),paste(mutations$st,mutations$wd,mutations$ty,sep="-")),file=paste(path,"/opLinks",path,".txt",sep=""),append=T,sep=",",col.names=F,row.names=F,quote=F);
      write.table(x=cbind(rep.int(x=contigNames[i],times=length(mutations$st)),as.character(mutations$st),as.character(mutations$wd),as.character(mutations$ty)),file=paste(path,"/opLinksii",path,".txt",sep=""),append=T,sep=",",col.names=F,row.names=F,quote=F);
    }
    write(x=start(subject(localAlign))-start(pattern(localAlign))+1,file=paste(path,"/startpositions",path,".txt",sep=""),ncolumns=1,append=T);  
  }
}

decomposeReadMatrix<-function(readFile,matrix,matrixFile,snpList,startPosition,contigLength,coverage,totalReads,path,contig,readLength)
{
  data <- matrix;
  reads <- read.delim(readFile, sep=',', header=FALSE, comment.char='>');
  
  subMatrixCount = 0;
  cumSumReads = 0;
  i = 1;
  
  while(i <= nrow(data))
  {
    j = 1;
    for(j in 1:(nrow(data)-(i-1)))
    {
      if((i+j) > nrow(data))
      {
        rightBlock = NULL;
      }
      else
      {
        rightBlock = data[i:(i+j-1),(i+j):ncol(data)];
      }  
      
      if(sum(rightBlock) == 0)
      {
        subMatrixCount = subMatrixCount + 1;
        subMatrix = data[i:(i+j-1),i:(i+j-1)];
        
        snpdec = matrix(snpList[i:(i+j-1),2:4],nrow=j,ncol=3,byrow=F);
        
        decName = paste(substr(matrixFile,1,nchar(matrixFile)-8),"_sub_",as.character(subMatrixCount),sep="");
        
        snpdec = cbind(decName,snpdec);
        
        diagonalSum = 0;  
        for(dim in 1:j){diagonalSum = diagonalSum + as.matrix(subMatrix)[dim,dim];}    
        sumReads = as.integer(0.5*(sum(subMatrix)+diagonalSum));
        
        if(sumReads > 0)
        {
          decStart = max(reads[cumSumReads+1,2],1) + startPosition - 1;
          cumSumReads = cumSumReads + sumReads;
          decEnd = max(reads[1:cumSumReads,3]) + startPosition - 1;         
          
          decLen = decEnd - decStart + 1;
          decCov = sumReads*readLength/decLen;
          details = c(decName,startPosition,contigLength,totalReads,coverage,decStart,decLen,sumReads,decCov);         
          
          decSeq = substr(x=contig,start=decStart-startPosition+1,stop=decEnd-startPosition+1);        
          
          write.table(x=snpdec,file=paste(path,"/opMutationsDec",path,".txt",sep=""),row.names=F,col.names=F,append=T,sep=",");
          write(x=details,file=paste(path,"/opclsters",path,".txt",sep=""),ncolumns=length(details),append=T,sep=",");
          write.fasta(sequences=decSeq,names=decName,file.out=paste(path,"/opclstersSeq",path,".fsa",sep=""),open="a",nbchar=decLen);      
          
          if((as.character(snpList[i,4])=="I")|(as.character(snpList[i,4])=="D"))
          {
            if(as.character(snpList[i,4])=="I")
            {
              stSnp = as.integer(snpList[i,2]);
            }
            if(as.character(snpList[i,4])=="D")
            {
              stSnp = as.integer(snpList[i,2])+1;
            }
          }
          else
          {
            stSnp = as.integer(snpList[i,2]);
          }
          
          if((as.character(snpList[(i+j-1),4])=="I")|(as.character(snpList[(i+j-1),4])=="D"))
          {
            if(as.character(snpList[(i+j-1),4])=="I")
            {
              enSnp = as.integer(snpList[(i+j-1),2])+1;
            }
            if(as.character(snpList[(i+j-1),4])=="D")
            {
              enSnp = as.integer(snpList[(i+j-1),2])+as.integer(nchar(snpList[(i+j-1),3]));
            }
          }
          else
          {
            enSnp = as.integer(snpList[(i+j-1),2]);
          }
          
          islandDetails = c(decName,stSnp,(enSnp-stSnp+1));
          
          write(x=islandDetails,file=paste(path,"/opIslands",path,".txt",sep=""),ncolumns=length(islandDetails),append=T,sep=",");
        }
        
        i = i + j;
        if(subMatrixCount > 100){return();}
        break;
      }
    }
  }
}

correctClusterError<-function(readFile,startPosition,contigLength,snpList,coverage,totalReads,path,contig,readLength,ctMutationPositions)
{
  data <- read.delim(readFile, sep=',', header=FALSE, comment.char='>');
  contigSP = startPosition;
  contigLn = contigLength;
  
  snps = ctMutationPositions[,2];
  snpmat = matrix(as.numeric(snps));
  
  readmat = matrix(data=0,nrow=nrow(snpmat),ncol=nrow(snpmat));
  for(i in 1:nrow(data))
  {
    for(j in 1:nrow(snpmat))
    {
      if (data[i,2]<=snpmat[j])
      {
        for (k in nrow(snpmat):j)
        {
          if (data[i,3]>=snpmat[k])
          {
            readmat[j,k]=readmat[j,k] + 1;
            readmat[k,j]=readmat[j,k];
            break;
          }
        }
        break;
      }
    } 
  }
  
  matrix=readmat;
  matrixFile = paste(substr(readFile,1,nchar(readFile)-4),"_mat.txt",sep="");
  decomposeReadMatrix(readFile,matrix,matrixFile,snpList,startPosition,contigLength,coverage,totalReads,path,contig,readLength);
}

estimateFreq<-function(path, clusterFile, alnStEndFile, startpos)
{
  clusterData = read.delim(file=clusterFile,header=F,sep=",");
  
  commonAlnData = read.delim(file=alnStEndFile,header=F,sep=",");
  
  filterAlnData = filterAln(path, startpos);
  
  write.table(x=filterAlnData,file=paste(path,"/filterReadStEnd",path,".txt",sep=""),row.names=F,col.names=F,sep=",");
  
  subNames = as.matrix(unlist(clusterData[,1]));
  subStartPos = as.matrix(as.integer(unlist(clusterData[,6])));
  subLength = as.matrix(as.integer(unlist(clusterData[,7])));
  subReadCount = as.matrix(as.integer(unlist(clusterData[,8])));
  
  commonReadCount = matrix(0,length(subNames),1);
  filterReadCount = matrix(0,length(subNames),1);
  freq = matrix(0,length(subNames),1)
  
  for (i in 1:length(subNames))
  {
    commonReadCount[i,1] = sum((commonAlnData[,1] >= subStartPos[i,1])&(commonAlnData[,2]<=subStartPos[i,1]+subLength[i,1]-1));
    filterReadCount[i,1] = sum((filterAlnData[,1] >= subStartPos[i,1])&(filterAlnData[,2]<=subStartPos[i,1]+subLength[i,1]-1));
    
  }
  freq = cbind(subNames,(subReadCount/(commonReadCount+filterReadCount)))
  
  for(i in 1:nrow(freq))
  {
    if(freq[i,2]>1)
      freq[i,2]=1;
  }
  
  write.table(x=freq,file=paste(path,"/opclstersFrq",path,".txt",sep=""),sep=",",append=F,row.names=F,col.names=F);
}

filterAln<-function(path, startpos)
{
  cts = c(1:nrow(startpos));
  filAlnStEnd = NULL;
  
  for(c in cts)
  {
    stpos = startpos[c];
    reads <- read.delim(paste(path,"/contig",as.character(c),".txt",sep=""), sep=',', header=FALSE, comment.char='>');
    alnpos = as.matrix(reads[,2]+stpos-1);
    alnend = as.matrix(reads[,3]+stpos-1);
    
    filAlnStEndNew = cbind(alnpos,alnend);
    filAlnStEnd = rbind(filAlnStEnd,filAlnStEndNew)
  } 
  
  return(filAlnStEnd);
}

overlapType<-function(s1,l1,s2,l2)
{ 
  if(((s1 <= s2) && (s2 <= (s1+l1-1))) || ((s2 <= s1) && (s1 <= (s2+l2-1))))
  {
    if((s1 == s2))
    {
      if(l1 == l2){return(1)}
      if(l1 <  l2){return(2)}
      if(l1 >  l2){return(3)}
    }
    if(s1 < s2)
    {
      if((s1+l1-1) == (s2+l2-1)){return(4)}
      if((s1+l1-1) <  (s2+l2-1)){return(5)}
      if((s1+l1-1) >  (s2+l2-1)){return(6)}
      
    }
    if(s1 > s2)
    {
      if((s1+l1-1) == (s2+l2-1)){return(7)}
      if((s1+l1-1) <  (s2+l2-1)){return(8)}
      if((s1+l1-1) >  (s2+l2-1)){return(9)}    
    }
  }
  else{return(0)}
}

offset<-function(ML,A,s,e,clsName)
{
  DL= ML[(as.character(A[,1])==clsName)&(as.integer(A[,2])>=s)&(as.integer(A[,2])<e)&(as.character(A[,4])=="D")]
  IL= ML[(as.character(A[,1])==clsName)&(as.integer(A[,2])>=s)&(as.integer(A[,2])<e)&(as.character(A[,4])=="I")]
  ofs = sum(DL)-sum(IL)
  return(ofs)
}

isOverlapAgreed<-function(str1,s1,l1,str2,s2,l2,clsName1,clsName2,mutationsData)
{ 
  offset1=0;
  offset2=0;
  A = mutationsData;
  ML = nchar(as.character(A[,3]));
  
  offset1 = sum(ML[as.character(A[as.character(A[,1])==clsName1,4])=="D"])-sum(ML[as.character(A[as.character(A[,1])==clsName1,4])=="I"])
  offset2 = sum(ML[as.character(A[as.character(A[,1])==clsName2,4])=="D"])-sum(ML[as.character(A[as.character(A[,1])==clsName2,4])=="I"])
  
  e1 = s1+l1+offset1-1
  e2 = s2+l2+offset2-1
  type = overlapType(s1,l1+offset1,s2,l2+offset2)
  
  if(type == 0)
  {
    return(2)
  }
  if(type == 1)
  {     
    if(str1 == str2){return(1)}
    else{return(0)}
  }
  if(type == 2)
  {
    poffset1s = offset(ML,A,s1,s1,clsName1)
    poffset1e = offset(ML,A,s1,e1,clsName1)
    
    poffset2s = offset(ML,A,s2,s1,clsName2)
    poffset2e = offset(ML,A,s2,e1,clsName2)
    
    ss1 = substr(x=str1,start=s1-poffset1s-s1+1,stop=e1-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s1-poffset2s-s2+1,stop=e1-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 3)
  {
    poffset1s = offset(ML,A,s1,s2,clsName1)
    poffset1e = offset(ML,A,s1,e2,clsName1)
    
    poffset2s = offset(ML,A,s2,s2,clsName2)
    poffset2e = offset(ML,A,s2,e2,clsName2)
    
    ss1 = substr(x=str1,start=s2-poffset1s-s1+1,stop=e2-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s2-poffset2s-s2+1,stop=e2-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 4)
  {
    poffset1s = offset(ML,A,s1,s2,clsName1)
    poffset1e = offset(ML,A,s1,e2,clsName1)
    
    poffset2s = offset(ML,A,s2,s2,clsName2)
    poffset2e = offset(ML,A,s2,e2,clsName2)
    
    ss1 = substr(x=str1,start=s2-poffset1s-s1+1,stop=e2-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s2-poffset2s-s2+1,stop=e2-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 5)
  {
    poffset1s = offset(ML,A,s1,s2,clsName1)
    poffset1e = offset(ML,A,s1,e1,clsName1)
    
    poffset2s = offset(ML,A,s2,s2,clsName2)
    poffset2e = offset(ML,A,s2,e1,clsName2)
    
    ss1 = substr(x=str1,start=s2-poffset1s-s1+1,stop=e1-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s2-poffset2s-s2+1,stop=e1-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 6)
  {
    poffset1s = offset(ML,A,s1,s2,clsName1)
    poffset1e = offset(ML,A,s1,e2,clsName1)
    
    poffset2s = offset(ML,A,s2,s2,clsName2)
    poffset2e = offset(ML,A,s2,e2,clsName2)
    
    ss1 = substr(x=str1,start=s2-poffset1s-s1+1,stop=e2-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s2-poffset2s-s2+1,stop=e2-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 7)
  {
    poffset1s = offset(ML,A,s1,s1,clsName1)
    poffset1e = offset(ML,A,s1,e1,clsName1)
    
    poffset2s = offset(ML,A,s2,s1,clsName2)
    poffset2e = offset(ML,A,s2,e1,clsName2)
    
    ss1 = substr(x=str1,start=s1-poffset1s-s1+1,stop=e1-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s1-poffset2s-s2+1,stop=e1-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 8)
  {
    poffset1s = offset(ML,A,s1,s1,clsName1)
    poffset1e = offset(ML,A,s1,e1,clsName1)
    
    poffset2s = offset(ML,A,s2,s1,clsName2)
    poffset2e = offset(ML,A,s2,e1,clsName2)
    
    ss1 = substr(x=str1,start=s1-poffset1s-s1+1,stop=e1-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s1-poffset2s-s2+1,stop=e1-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
  if(type == 9)
  {
    poffset1s = offset(ML,A,s1,s1,clsName1)
    poffset1e = offset(ML,A,s1,e2,clsName1)
    
    poffset2s = offset(ML,A,s2,s1,clsName2)
    poffset2e = offset(ML,A,s2,e2,clsName2)
    
    ss1 = substr(x=str1,start=s1-poffset1s-s1+1,stop=e2-poffset1e-s1+1);
    ss2 = substr(x=str2,start=s1-poffset2s-s2+1,stop=e2-poffset2e-s2+1);
    if(ss1 == ss2){return(1)}
    else{return(0)}
  }
}

createMutationGraph<-function(path,clusterFile,islandsFile,freqFile,stringFile,refLength)
{  
  freqFile = paste(path,"/opclstersFrq",path,".txt",sep="");
  freqData = read.delim(file=freqFile,header=F,sep=",");
  fN = as.character(freqData[,1]);
  fF = as.numeric(freqData[,2]);
  
  islandsFile=paste(path,"/opIslands",path,".txt",sep="");
  islandsData = read.delim(file=islandsFile,header=F,sep=",");
  subNames = as.matrix((unlist(islandsData[,1])));
  subStartPos = as.matrix(as.integer(unlist(islandsData[,2])));
  subLength = as.matrix(as.integer(unlist(islandsData[,3])));
  
  mutationsFile = paste(path,"/opMutationsDec",path,".txt",sep="");
  mutationsData = read.delim(file=mutationsFile,header=F,sep=",");
  
  clusterFile = paste(path,"/opclsters",path,".txt",sep="");
  clusterData = read.delim(file=clusterFile,header=F,sep=",");
  clsNames = as.matrix((unlist(clusterData[,1])));
  clsStartPos = as.matrix(as.integer(unlist(clusterData[,6])));
  clsLength = as.matrix(as.integer(unlist(clusterData[,7])));
  clsRNs = as.matrix(as.integer(unlist(clusterData[,8])));
  
  referenceString = unlist(read.fasta(file="SIM/reference.fsa",as.string=T,seqonly=T));
  stringFile = paste(path,"/opclstersSeq",path,".fsa",sep="");
  strings = unlist(read.fasta(file=stringFile,as.string=T,seqonly=T));
  
  overlapMat = matrix(data=0,nrow=length(strings),ncol=length(strings));
  mutOverlapMat = matrix(data=0,nrow=length(strings),ncol=length(strings));
  toberem = c();
  extended = c();
  
  for(i in 1:length(strings))
  {
    for(j in i:length(strings))
    {
      overlapMat[i,j]=(isOverlapAgreed(strings[i],clsStartPos[i,1],clsLength[i,1],
                                       strings[j],clsStartPos[j,1],clsLength[j,1],
                                       clsNames[i,1],clsNames[j,1],mutationsData) == 1)
      overlapMat[j,i]=overlapMat[i,j]
      a = unique(paste(mutationsData[which(mutationsData[,1] == clsNames[i,1]),2],
                       mutationsData[which(mutationsData[,1] == clsNames[i,1]),3],
                       mutationsData[which(mutationsData[,1] == clsNames[i,1]),4],sep="-"));
      b = unique(paste(mutationsData[which(mutationsData[,1] == clsNames[j,1]),2],
                       mutationsData[which(mutationsData[,1] == clsNames[j,1]),3],
                       mutationsData[which(mutationsData[,1] == clsNames[j,1]),4],sep="-"));
      c = unique(c(a,b));
      mutOverlapMat[i,j]=((length(a)+length(b))>length(c))
      mutOverlapMat[j,i]=mutOverlapMat[i,j]     
      
      if(i != j)
      {
        if(overlapMat[i,j]&mutOverlapMat[i,j])
        {
          
          snpList = mutationsData[(as.character(mutationsData[,1])==clsNames[i,1])|(as.character(mutationsData[,1])==clsNames[j,1]),]
          snpList = snpList[order(snpList[,2]),]
          if((as.character(snpList[1,4])=="I")|(as.character(snpList[1,4])=="D"))
          {
            if(as.character(snpList[1,4])=="I")
            {
              stSnp = as.integer(snpList[1,2]);
            }
            if(as.character(snpList[1,4])=="D")
            {
              stSnp = as.integer(snpList[1,2])+1;
            }
          }
          else
          {
            stSnp = as.integer(snpList[1,2]);
          }
          
          if((as.character(snpList[nrow(snpList),4])=="I")|(as.character(snpList[nrow(snpList),4])=="D"))
          {
            if(as.character(snpList[nrow(snpList),4])=="I")
            {
              enSnp = as.integer(snpList[nrow(snpList),2])+1;
            }
            if(as.character(snpList[nrow(snpList),4])=="D")
            {
              enSnp = as.integer(snpList[nrow(snpList),2])+as.integer(nchar(as.character(snpList[nrow(snpList),3])));
            }
          }
          else
          {
            enSnp = as.integer(snpList[nrow(snpList),2]);
          }
          
          subStartPos[i,1] = stSnp;
          subLength[i,1] = enSnp-stSnp+1;
          
          nms = c(as.character(clusterData[i,1]),as.character(clusterData[j,1]))
          strs = c(strings[i],strings[j])
          starts = c(clusterData[i,6],clusterData[j,6]);
          lengths = c(clusterData[i,7],clusterData[j,7]);
          ends = c((clusterData[i,6]+clusterData[i,7]-1),(clusterData[j,6]+clusterData[j,7]-1));
          st = min(starts);
          en = max(ends);
          
          startcon = which(starts==st)[1];
          endcon = which(ends==en)[1];
          
          A = mutationsData;
          ML = nchar(as.character(A[,3]));
          offset1 = sum(ML[as.character(A[as.character(A[,1])==nms[startcon],4])=="D"])-sum(ML[as.character(A[as.character(A[,1])==nms[startcon],4])=="I"])
          offset2 = sum(ML[as.character(A[as.character(A[,1])==nms[endcon],4])=="D"])-sum(ML[as.character(A[as.character(A[,1])==nms[endcon],4])=="I"])
          
          e1 = starts[startcon]+lengths[startcon]+offset1-1
          e2 = starts[endcon]+lengths[endcon]+offset2-1
          
          poffset1e = offset(ML,A,starts[startcon],e1,nms[startcon])
          poffset2e = offset(ML,A,starts[endcon],e1,nms[endcon])
          
          strings[i] = paste(substr(strs[startcon],start=1,stop=lengths[startcon]),
                             substr(strs[endcon],start=starts[startcon]+lengths[startcon]+poffset1e-poffset2e-starts[endcon]+1,stop=lengths[endcon]),sep="");
          
          clsStartPos[i,1] = st;
          clsLength[i,1] = en-st+1;
          clsRNs[i,1] = clsRNs[i,1]+clsRNs[j,1]
          
          mutationsData[which(mutationsData[,1] == clsNames[j,1]),1] = clsNames[i,1];
          
          toberem = c(toberem,j)
          extended = c(extended,i)
        }
      }
    }
  }
  
  write.table(x=cbind(0,0),file=paste(path,"/combined.txt",sep=""),append=F,sep="\t",row.names=F,col.names=F)
  
  
  if (length(toberem)>0)
  {
    write.table(x=cbind(clsNames[extended,1],clsNames[toberem,1]),file=paste(path,"/combined.txt",sep=""),append=F,sep="\t",row.names=F,col.names=F)
    
    rowNum = 1:nrow(islandsData);
    toberem = sort(unique(toberem));
    for(i in toberem)
    {
      rowNum = rowNum[-which(rowNum == i)]
    }
    
    subNames = subNames[rowNum,1];
    subStartPos = subStartPos[rowNum,1];
    subLength = subLength[rowNum,1];
    
    clsNames = clsNames[rowNum,1];
    clsStartPos = clsStartPos[rowNum,1];
    clsLength = clsLength[rowNum,1];
    clsRNs = clsRNs[rowNum,1];
    
    strings = strings[rowNum];
  }  
  
  commonAlnData = read.delim(file=paste(path,"/commonReadAlnStEnd",path,".txt",sep=""),header=F,sep=",");  
  filterAlnData = read.delim(file=paste(path,"/filterReadStEnd",path,".txt",sep=""),header=F,sep=","); 
  commonReadCount = matrix(0,length(clsNames),1);
  filterReadCount = matrix(0,length(clsNames),1);
  freq = matrix(0,length(clsNames),1)
  
  for (i in 1:length(clsNames))
  {
    commonReadCount[i,1] = sum((as.integer(commonAlnData[,1]) >= clsStartPos[i])&(as.integer(commonAlnData[,2])<=clsStartPos[i]+clsLength[i]-1));
    filterReadCount[i,1] = sum((as.integer(filterAlnData[,1]) >= clsStartPos[i])&(as.integer(filterAlnData[,2])<=clsStartPos[i]+clsLength[i]-1));
    
  }
  freq = cbind(clsNames,(clsRNs/(commonReadCount+filterReadCount)))
  fN = as.character(freq[,1]);
  for(i in 1:nrow(freq))
  {
    if(freq[i,2]>1)
      freq[i,2]=1;
  }
  fF = as.numeric(freq[,2]);
  
  newislandsFile = paste(path,"/newopIslands",path,".txt",sep="");
  write.table(x=cbind(subNames,subStartPos,subLength),file=newislandsFile,append=F,row.names=F,col.names=F,sep=",")
  
  newclusterFile = paste(path,"/newopclsters",path,".txt",sep="");
  write.table(x=cbind(clsNames,clsStartPos,clsLength,clsRNs),file=newclusterFile,append=F,row.names=F,col.names=F,sep=",")
  
  newmutationsFile = paste(path,"/newopMutationsDec",path,".txt",sep="");
  write.table(x=mutationsData,file=newmutationsFile,append=F,row.names=F,col.names=F,sep=",")
  
  sMat = cbind(subNames,"s",subStartPos-0.5);
  eMat = cbind(subNames,"e",subStartPos+subLength-1+0.5);
  
  seMat = rbind(sMat,eMat);
  
  seMat = seMat[order(as.integer(seMat[,3])),]
  
  mN = seMat[,1]; mSE = seMat[,2]; mL = as.numeric(seMat[,3]);
  seMat = data.frame(mN,mSE,mL);
  
  unmL = unique(seMat$mL)
  
  nMut = nrow(sMat);
  
  nRef = length(unmL) + 1;
  
  for(i in 0:(nRef - 1))
  {
    if(i == 0)
    {
      mN = c(mN,paste(path,"/","ref",as.character(i),sep=""),paste(path,"/","ref",as.character(i),sep=""));
      mSE = c(mSE,"s","e");
      mL = c(mL,1,unmL[i+1]);
      fN = c(fN,paste(path,"/","ref",as.character(i),sep=""));
      fF = c(fF,0)
      
      clsNames = c(clsNames,paste(path,"/","ref",as.character(i),sep=""))
      clsStartPos = c(clsStartPos,1)
      clsLength = c(clsLength,((unmL[i+1]-0.5)-1+1))
      clsRNs = c(clsRNs,0)
      
      refstr = substr(x=referenceString,start=1,stop=(unmL[i+1]-0.5))
      strings = c(strings,refstr)
    }
    if(i == length(unmL))
    {
      mN = c(mN,paste(path,"/","ref",as.character(i),sep=""),paste(path,"/","ref",as.character(i),sep=""));
      mSE = c(mSE,"s","e");
      mL = c(mL,unmL[i],refLength);
      fN = c(fN,paste(path,"/","ref",as.character(i),sep=""));
      fF = c(fF,0)
      
      clsNames = c(clsNames,paste(path,"/","ref",as.character(i),sep=""))
      clsStartPos = c(clsStartPos,unmL[i]+0.5)
      clsLength = c(clsLength,(refLength-(unmL[i]+0.5)+1))
      clsRNs = c(clsRNs,0)
      
      refstr = substr(x=referenceString,start=(unmL[i]+0.5),stop=refLength)
      strings = c(strings,refstr)
    }
    if((i > 0) && (i < length(unmL)))
    {
      mN = c(mN,paste(path,"/","ref",as.character(i),sep=""),paste(path,"/","ref",as.character(i),sep=""));
      mSE = c(mSE,"s","e");
      mL = c(mL,unmL[i],unmL[i+1]);
      fN = c(fN,paste(path,"/","ref",as.character(i),sep=""));
      fF = c(fF,0);
      
      clsNames = c(clsNames,paste(path,"/","ref",as.character(i),sep=""))
      clsStartPos = c(clsStartPos,unmL[i]+0.5)
      clsLength = c(clsLength,((unmL[i+1]-0.5)-(unmL[i]+0.5)+1))
      clsRNs = c(clsRNs,0)
      
      refstr = substr(x=referenceString,start=(unmL[i]+0.5),stop=(unmL[i+1]-0.5))
      strings = c(strings,refstr)
    }
  }
  #sourse
  mN = c(mN,paste(path,"/S",sep=""),paste(path,"/S",sep=""));
  mSE = c(mSE,"s","e");
  mL = c(mL,0.5,1);
  #sink
  mN = c(mN,paste(path,"/T",sep=""),paste(path,"/T",sep=""));
  mSE = c(mSE,"s","e");
  mL = c(mL,refLength,refLength+0.5);
  
  seMat = data.frame(mN,mSE,mL);
  seMat = seMat[order(as.integer(seMat[,3]),seMat[,2]),];
  
  tailMat = seMat[(seMat[,2]=="e"),]
  headMat = seMat[(seMat[,2]=="s"),]
  
  tailMat = tailMat[order(tailMat$mN),]
  headMat = headMat[order(headMat$mN),]
  
  freq = data.frame(fN,fF);
  freq = freq[order(freq$fN),];
  
  nodeOrder = as.character(headMat$mN);
  
  adjacencyMat = matrix(data=0,nrow=nrow(tailMat),ncol=nrow(headMat),dimnames=list(c(as.character(tailMat$mN)),c(as.character(headMat$mN))));
  adjacencyList = matrix(data=c(1:(nrow(tailMat))),nrow=(nrow(tailMat)),ncol=2)
  
  for(i in 1:nrow(tailMat))
  {
    adjacencyMat[i,] = as.integer(headMat[,3] == tailMat[i,3]);
    adjacencyList[i,2] = paste((which(headMat[,3] == tailMat[i,3])),collapse=",");
  }
  
  adjacencyListFile = paste(path,"/adjacencyList",path,".txt",sep="");  
  write.table(x=adjacencyList,file=adjacencyListFile,append=F,sep="-",row.names=F,col.names=F,quote=F);
  
  adjacencyFile = paste(path,"/adjacencyMatrix",path,".txt",sep="");
  write.table(x=adjacencyMat,file=adjacencyFile,append=F,sep=",");
  
  constMat = matrix(data=0,nrow=nrow(tailMat)-2,ncol=nrow(headMat)-2,dimnames=list(c(1:(length(tailMat$mN)-2)),c(1:(length(headMat$mN)-2))));
  overlapNumMat = matrix(data=0,nrow=nrow(tailMat)-2,ncol=1,dimnames=list(c(1:(length(tailMat$mN)-2)),"OverlapNum"))
  
  for(i in 1:(nrow(tailMat)-2))
  {
    constMat[i,] = 1 - as.integer(((headMat[,3] >= (tailMat[i,3]+0))|(tailMat[,3] <= (headMat[i,3]-0)))[1:(nrow(headMat)-2)]);
    overlapNumMat[i,1] = sum(constMat[i,]);
  }    
  
  write.table(x=headMat[1:(nrow(headMat)-2),],file=paste(path,"/headMat.txt",sep=""),sep=",",quote=F,row.names=F,col.names=F);
  write.table(x=clsLength,file=paste(path,"/clsLenMat.txt",sep=""),sep=",",quote=F,row.names=F,col.names=F);
  write.table(x=constMat,file=paste(path,"/constMat.txt",sep=""),sep=",",quote=F,row.names=F,col.names=F);
  write.table(x=nRef,file=paste(path,"/nRef.txt",sep=""),sep=",",quote=F,row.names=F,col.names=F);
  
  overlapMat = matrix(data=0,nrow=nrow(tailMat)-2,ncol=nrow(headMat)-2,dimnames=list(c(1:(length(tailMat$mN)-2)),c(1:(length(headMat$mN)-2))));
  
  newclusterData = data.frame(clsNames,clsStartPos,clsLength,clsRNs,strings)
  newclusterData = newclusterData[order(newclusterData$clsNames),]
  
  for(i in 1:nrow(newclusterData))
  {
    for(j in i:nrow(newclusterData))
    {
      overlapMat[i,j]=isOverlapAgreed(as.character(newclusterData[i,5]),newclusterData[i,2],newclusterData[i,3],
                                      as.character(newclusterData[j,5]),newclusterData[j,2],newclusterData[j,3],
                                      as.character(newclusterData[i,1]),as.character(newclusterData[j,1]),mutationsData)
      overlapMat[j,i]=overlapMat[i,j]
      
    }
  }
  overlapMatrixFile = paste(path,"/overlapMatrix",path,".txt",sep="");
  write.table(x=overlapMat,file=overlapMatrixFile,append=F,sep=",",row.names=F,col.names=F);
  
  for(i in (nrow(freq)-nRef+1):nrow(freq))
  {
    freq[i,2] = min((1 - constMat[i,] %*% freq$fF),1);
  }
  
  freqListFile = paste(path,"/freqList",path,".txt",sep="");
  
  zfreqs = which(freq[,2]<0);
  if(length(zfreqs)>0)
  {
    freq[zfreqs,2] = 0;
  }  
  write.table(x=freq[,2],file=freqListFile,append=F,sep=",",row.names=F,col.names=F);
}

constructStrainsDynamic<-function(path,mode,cutoff){
  
  bp="0";
  
  ############# Initialize #########################################################
  
  headMat   = read.table(file=paste(path,"/headMat.txt",sep=""),header=F,sep=",",quote="");
  lengthMat = read.table(file=paste(path,"/clsLenMat.txt",sep=""),header=F,sep=",",quote="");
  constMat  = read.table(file=paste(path,"/constMat.txt",sep=""),header=F,sep=",",quote="",stringsAsFactors=F);
  freqList  = read.table(file=paste(path,"/freqList",path,".txt",sep=""),header=F,sep=",",quote="");
  nRef      = read.table(file=paste(path,"/nRef.txt",sep=""),header=F,sep=",",quote="")[1,1];
  overlapMat= read.table(file=paste(path,"/overlapMatrix",path,".txt",sep=""),header=F,sep=",",quote="");
  
  partitionOrder = order(headMat[(nrow(headMat)-nRef+1):nrow(headMat),3]);
  partitionNames = headMat[(nrow(headMat)-nRef+1):nrow(headMat),1];
  partitionNames = partitionNames[partitionOrder];
  
  subNames  = NULL;subName   = NULL;subFreq   = NULL;subPart   = NULL;subLength= NULL;
  
  for(c in (nrow(headMat)-nRef+1):nrow(headMat)){
    newsubNames= which(constMat[c,]==1);
    subName    = c(subName  ,newsubNames)
    subNames   = c(subNames ,as.character(headMat[newsubNames,1]))
    subFreq    = c(subFreq  ,freqList[newsubNames,1]/sum(freqList[newsubNames,1]))
    #subFreq    = c(subFreq  ,freqList[newsubNames,1])
    subPart    = c(subPart  ,rep.int(x=which(partitionNames==headMat[c,1]),times=length(newsubNames)))
    subLength  = c(subLength,lengthMat[newsubNames,1])
  }
  
  partitionData = data.frame(subName,subFreq,subPart,subNames,subLength,subFreq,stringsAsFactors=F);
  
  nPartitions      = nRef;
  constructedPaths = NULL;
  constructedFreqs = NULL;
  newIdx           = matrix(data=0,nrow=1,ncol=nPartitions);
  
  maxData       = matrix(data=0,nrow=nrow(constMat),ncol=1);
  for(i in 1:nrow(constMat)){maxData[i,1] = sum(constMat[i,])}
  maxIterations = 2000;
  
  allZeros = FALSE; i = 1;
  
  ############# Iterations #########################################################
  
  while((!allZeros)&(i <= maxIterations))
  {   
    guide = guidePartition(partitionData);
    
    newPath = matrix(data=0,nrow=2,ncol=nPartitions);
    m1 = (partitionData[,3]==guide);
    max = max(partitionData[m1,2]);
    m2 = (partitionData[,2]==max);
    m3 = m1&m2;
    
    newIdx[1,] = which(m3)[1];
    if(is.na(newIdx[1,guide])){bp="1";break;}
    newPath[1,guide] = as.integer(partitionData[newIdx[1,guide],1]);
    newPath[2,guide] = partitionData[newIdx[1,guide],2];
    
    if(guide < ncol(newPath))
    {     
      for(j in (guide+1):ncol(newPath))
      {
        
        m1 = (partitionData[,3]==j); 
        m2 = (partitionData[,1]==newPath[1,j-1]);
        m3 = m1&m2;
        
        common = (length(which(m3))>0);
        if(common)
        {
          newIdx[1,j] = which(m3)[1];
          if(is.na(newIdx[1,j])){bp="2";break;}
          newPath[1,j] = as.integer(partitionData[newIdx[1,j],1]);      
          newPath[2,j] = partitionData[newIdx[1,j],2];
          
        }
        else
        {
          distance = abs(partitionData[,2]-as.numeric(min(newPath[2,newPath[2,]>0])));
          
          overlapAgreed = c(1:nrow(partitionData));
          constrained   = c(1:nrow(partitionData));
          
          a= overlapMat[as.integer(newPath[1,guide:(j-1)]),partitionData[,1]]==0
          a = !a
          for(l in 1:nrow(partitionData))
          {
            overlapAgreed[l] = (sum(a[,l])>0)
            constrained[l]   = constMat[as.integer(newPath[1,j-1]),partitionData[l,1]]
          }
          
          oc = (overlapAgreed)&(!constrained);
          
          if(sum(m1&oc)<1){bp="3";break;}
          
          minDistance = min(distance[m1&oc]);
          m2 = (distance==minDistance);
          m3 = (m1&m2&oc);
          
          newIdx[1,j] = which(m3)[1];
          if(is.na(newIdx[1,j])){bp="4";break;}
          newPath[1,j] = as.integer(partitionData[newIdx[1,j],1]);        
          newPath[2,j] = partitionData[newIdx[1,j],2];
        }   
      }
    }
    
    if(guide > 1)
    {
      for(j in (guide-1):1)
      { 
        m1 = (partitionData[,3]==j); 
        m2 = (partitionData[,1]==newPath[1,j+1]);
        m3 = m1&m2;
        
        common = (length(which(m3))>0);
        if(common)
        {
          newIdx[1,j] = which(m3)[1];
          if(is.na(newIdx[1,j])){bp="5";break;}
          newPath[1,j] = as.integer(partitionData[newIdx[1,j],1]);
          newPath[2,j] = partitionData[newIdx[1,j],2];          
        }
        else
        {
          distance = abs(partitionData[,2]-as.numeric(min(newPath[2,newPath[2,]>0])));
          
          overlapAgreed = c(1:nrow(partitionData));
          constrained   = c(1:nrow(partitionData));
          
          a= overlapMat[as.integer(newPath[1,(j+1):guide]),partitionData[,1]]==0
          a = !a
          for(l in 1:nrow(partitionData))
          {
            overlapAgreed[l] = (sum(a[,l])>0)
            constrained[l]   = constMat[as.integer(newPath[1,j+1]),partitionData[l,1]]
          }
          
          oc = (overlapAgreed)&(!constrained);
          
          if(sum(m1&oc)<1){bp="6";break;}
          
          minDistance = min(distance[m1&oc]);
          m2 = (distance==minDistance);
          m3 = (m1&m2&oc);
          
          newIdx[1,j] = which(m3)[1];
          if(is.na(newIdx[1,j])){bp="7";break;}
          newPath[1,j] = as.integer(partitionData[newIdx[1,j],1]);         
          newPath[2,j] = partitionData[newIdx[1,j],2];
        }   
      }
    }
    
    if(!allZeros)
    {
      if(i==1)
      {
        constructedPaths = newPath[1,];
        constructedFreqs = newPath[2,];
      }
      else
      {
        constructedPaths = rbind(constructedPaths,newPath[1,]);
        constructedFreqs = rbind(constructedFreqs,newPath[2,]);
      }
      i = i + 1;
      
      if(mode == 1)
      {
        partitionData[newIdx,2] = as.numeric(partitionData[newIdx,2]) - as.numeric(min(newPath[2,(newPath[2,]>0)]));
      }
      if(mode == 2)
      {
        partitionData[newIdx,2] = as.numeric(partitionData[newIdx,2]) - as.numeric(partitionData[newIdx[guide],2]);       
      }
      
      
      if(sum(is.na(partitionData[,2]))>0)
      {
        bp = "8"
        break;
      }  
      
      partitionData[partitionData[,2]<0,2] = 0;
      
      for(k in 1:nPartitions)
      {
        if(sum(partitionData[which(partitionData[,3]==k),2])>0)
        {
          allZeros = allZeros|F;
        }
        else
        {
          bp = "9"
          allZeros = T;
        }
      }
      
      partitionData=partitionData[(partitionData[,2] > 0),];
      
    }
  }  
  
  write(x=t(constructedPaths),file=paste(path,"/constructedPaths",path,".txt",sep=""),append=F,ncolumns=nRef);
  write(x=t(constructedFreqs),file=paste(path,"/constructedFreqs",path,".txt",sep=""),append=F,ncolumns=nRef);
  newmutationsFile = paste(path,"/newopMutationsDec",path,".txt",sep="");
  opMutMat = read.table(file=newmutationsFile,header=F,sep=",");
  
  df <- try(constructedFreqs <- matrix(data=constructedFreqs,nrow=as.integer(length(constructedFreqs)/nPartitions),ncol=nPartitions,byrow=F),silent=TRUE)
  if (class(df)=='try-error'){  constructedFreqs = matrix(1,1,1) }
  
  fPaths = matrix(data=0,nrow=nrow(constructedFreqs),ncol=1);
  for(i in 1:nrow(constructedFreqs)){
    fPaths[i,1] = max(0,min(constructedFreqs[i,constructedFreqs[i,]>0]))
  }
  
  if(sum(fPaths)>1)
  {
    fPaths[,1] = fPaths[,1]/sum(fPaths[,1]);
  } 
  
  constructSequence(path);
  
  s = unlist(read.fasta(file=paste(path,"/reconstructedSeqs",path,".fsa",sep=''),as.string=T,forceDNAtolower=F,seqonly=T))
  b = unique(s);
  f = fPaths[!duplicated(s),1]
  
  for(i in 1:length(f))
  {
    f[i] = sum(fPaths[s==b[i],1])
  }
  
  if((length(b)==1)&(b[1]=="")){b="0"}
  
  resTable = cbind(f,b)
  if(length(f)>1)
  {
    resTable = resTable[order(as.numeric(resTable[,1]),decreasing=T),]  
  }
  
  f = as.numeric(resTable[,1]);
  b = resTable[,2]
  
  write(x=f,file=paste(path,"/reconstructedFreqs",path,".txt",sep=""),ncolumns=1,append=F)
  
  for(i in 1:length(b))
  {
    write.fasta(sequences=b[i],names=paste("Strain",as.character(i),"_",as.character(f[i]),sep=""),file.out="ViQuaS-Spectrum.fa",open="a",nbchar=60)   
  }
}

guidePartition<-function(partitionData){
  
  if(is.null(partitionData))
  {
    return(0);
  }
  else
  {
    partitionHist = hist(partitionData[,3],breaks=(0.5:(max(partitionData[,3])+0.5)),plot=F);
    guide = which(partitionHist$counts == max(partitionHist$counts))[1];
    return(guide);
  } 
}

constructSequence<-function(path){
  
  reference = read.fasta(file="SIM/reference.fsa", seqtype="DNA",forceDNAtolower=F,as.string=T);
  refString = (reference$reference)[1];
  
  headMat   = read.table(file=paste(path,"/headMat.txt",sep=""),header=F,sep=",",quote="");
  constructedPaths = read.table(file=paste(path,"/constructedPaths",path,".txt",sep=""),header=F,sep=" ");
  newmutationsFile = paste(path,"/newopMutationsDec",path,".txt",sep="");
  opMutMat = read.table(file=newmutationsFile,header=F,sep=",");
  
  constructedSequences = matrix(data="",nrow=nrow(constructedPaths),ncol=1)
  
  for(i in 1:nrow(constructedPaths)){
    refArray = unlist(strsplit(refString, split= ''));  
    nodelist = unique(as.character(headMat[unlist(constructedPaths[i,]),1]));
    mutationList = c();
    for (j in 1:length(nodelist))
    {
      mutationList = rbind(mutationList,opMutMat[which(as.character(opMutMat[,1])==nodelist[j]),2:4]);
    }   
    
    snpList = mutationList[(mutationList[,3]!="I")&(mutationList[,3]!="D"),]
    if(nrow(snpList)>1)
    {snpList = snpList[!duplicated(snpList[,1]),]}
    
    delList = mutationList[mutationList[,3]=="D",]
    if(nrow(delList)>1)
    {delList = delList[!duplicated(delList[,1]),]}
    
    insList = mutationList[mutationList[,3]=="I",]
    if(nrow(insList)>1)
    {insList = insList[!duplicated(insList[,1]),]}
    
    if(nrow(snpList) > 0)
    {
      for(j in 1:nrow(snpList)){refArray[snpList[j,1]] = as.character(snpList[j,3]);}
    } 
    if(nrow(delList) > 0)
    {
      for(j in 1:nrow(delList))
      {
        for(k in 1:nchar(as.character(delList[j,2])))
        {
          refArray[delList[j,1]+k] = NA;
        }
      }
    }
    if(nrow(insList) > 0)
    {
      for(l in nrow(insList):1)
      {
        refArray = c(refArray[1: (as.integer(insList[l,1])-0)],unlist(strsplit(as.character(insList[l,2]),'')),refArray[(as.integer(insList[l,1])+1): as.integer(length(refArray))]);
      }
    }    
    write.fasta(sequences=na.exclude(refArray),names=as.character(i),file.out=paste(path,"/reconstructedSeqs",path,".fsa",sep=''),open="a");
  }
}

################# QNE ################################################################################

QNE<-function(path,refLength,readLength,nReads,diversityRegionLength)
{
  cutoff = (-1)*refLength*log(1-0.99^(1/refLength))/readLength/nReads;
  #cutoff = round(cutoff,4);
  
  clFreqs = read.table(file=paste(path,"/opclstersFrq",path,".txt",sep=""),header=F,sep=',',stringsAsFactors=F);
  clf = clFreqs[clFreqs[,2]>cutoff,1]
  clm = read.table(file=paste(path,"/newopMutationsDec",path,".txt",sep=""),header=F,sep=',',stringsAsFactors=T);
  opm = unique(apply(subset(x=clm,subset= clm$V1 %in% clf)[,2:4],1,paste,collapse="-"))
  
  l = strsplit(x=opm,split="-");
  m=c();
  for(u in 1:length(l))
  {
    m = c(m,nchar(l[[u]][2]))
  }
  opnuMutp = sum(as.numeric(m))
  
  ############################################################
  
  reads = readDNAStringSet(paste("Samples/",path,"/filtered.fsa",sep=""), "fasta");
  reference = readDNAStringSet("SIM/reference.fsa", "fasta");
  numReads = length(reads);
  
  allmutations = c();
  
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE);
  
  for(a in 1:numReads)
  {
    localAlign = pairwiseAlignment(reads[a], reference[1], type="overlap" ,substitutionMatrix = mat, gapOpening = -5, gapExtension = -2);
    
    snps = mismatchSummary(localAlign);
    st = as.integer(snps$subject$SubjectPosition);
    wd = c(rep(1,length(snps$subject$SubjectPosition)));
    ty = snps$subject$Pattern;
    opsnps = data.frame(st,wd,ty);
    st = as.integer(snps$subject$SubjectPosition) - start(subject(localAlign)) + start(pattern(localAlign));
    opsnpsCon = data.frame(st,wd,ty);
    
    loc = NULL;
    locCon = NULL;
    locRef = NULL;
    
    indels = indel(localAlign);
    ins = indels@insertion@unlistData;     
    del = indels@deletion@unlistData;
    
    if(length(ins)+length(del) > 0)
    {
      insRows = NULL;
      delRows = NULL;
      
      if(length(ins) > 0)
      {
        insRows = cbind(ins@start,ins@width, "I")
      }
      
      if(length(del) > 0)
      {
        delRows = cbind(del@start,del@width, "D")
      }
      
      locs = rbind(insRows,delRows)
      st = as.integer(locs[,1]);wd = as.integer(locs[,2]);ty = locs[,3];sq="";
      loc = data.frame(st,wd,ty,as.character(sq),stringsAsFactors=F)
      loc = loc[order(loc$st),]
      locRef = data.frame(st,wd,ty,as.character(sq),stringsAsFactors=F)
      locRef = locRef[order(locRef$st),]
      locCon = data.frame(st,wd,ty,as.character(sq),stringsAsFactors=F)
      locCon = locCon[order(locCon$st),]
      insoffset = 0;
      deloffset = 0;
      
      for(j in 1:nrow(locs))
      {
        if(j>1)
        {
          if(loc[j-1,3] == "I")
          {
            insoffset = insoffset + as.integer(loc[j-1,2]);
          }         
          if(locRef[j-1,3] == "D")
          {
            deloffset = deloffset + as.integer(loc[j-1,2]);
          }
        }         
        if(locRef[j,3] == "D")
        {
          locRef[j,1] = as.integer(locRef[j,1]) + start(subject(localAlign)) - start(pattern(localAlign)) - insoffset + deloffset - 1;
          locCon[j,1] = as.integer(locCon[j,1]) - 1;  
          
          mu = substr(x=reference[1],start=locRef[j,1]+1,stop=locRef[j,1]+as.integer(locRef[j,2]))
          locRef[j,4] = as.character(mu)
          locCon[j,4] = as.character(mu)
        } 
        if(locRef[j,3] == "I")
        {
          locRef[j,1] = as.integer(locRef[j,1]) + start(subject(localAlign)) - start(pattern(localAlign)) - 1;
          locCon[j,1] = as.integer(locCon[j,1]) + insoffset - deloffset - 1;
          
          mu = substr(x=reads[a],start=locCon[j,1]+1,stop=locCon[j,1]+as.integer(locCon[j,2]))
          locRef[j,4] = as.character(mu)
          locCon[j,4] = as.character(mu)
        }
      }
    }
    
    st = locCon$st;wd = locCon$as.character.sq.;ty=locCon$ty;
    contigPositions = rbind(opsnpsCon,cbind(st,wd,ty));
    contigPositions = contigPositions[order(as.integer(contigPositions$st)),];
    
    st = locRef$st;wd = locRef$as.character.sq.;ty=locRef$ty;
    mutations = rbind(opsnps,cbind(st,wd,ty));
    mutations = mutations[order(as.integer(mutations$st)),]; 
    
    allmutations = rbind(allmutations,mutations)
  }
  
  mutlocs = apply(X=allmutations,MARGIN=1,FUN=paste,collapse="-")
  
  l = strsplit(x=mutlocs,split="-");
  m=c();
  
  for(u in 1:length(l))
  {
    m = c(m,as.numeric(l[[u]][1]))
  }
  
  #totalMut = sum(as.numeric(m))
  
  totalMut = length(m)
  
  #   ml = unique(mutlocs)
  #   l = strsplit(x=ml,split="-")
  #   m=c();
  #   for(u in 1:length(l))
  #   {
  #     m = c(m,as.numeric(l[[u]][1]))
  #   }
  #opnuMut = sum(as.numeric(m))
  
  opnuMut = length(unique(m))
  
  ####################################################################################################
  
  strainsAll = readDNAStringSet(paste("ViQuaS-Spectrum.fa",sep=""), "fasta");  
  stFreqs = read.table(file=paste(path,"/reconstructedFreqs",path,".txt",sep=""),,header=F,sep=',',stringsAsFactors=F);
  strains = strainsAll[which(stFreqs[,1]>cutoff)]
  
  reference = readDNAStringSet("SIM/reference.fsa", "fasta");
  numStrains = length(strains);
  
  mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3, baseOnly = TRUE);
  
  nMuts = c();
  
  for(a in 1:numStrains)
  {
    localAlign = pairwiseAlignment(strains[a], reference[1], type="overlap" ,substitutionMatrix = mat, gapOpening = -5, gapExtension = -2);
    snps = mismatchSummary(localAlign);
    indels = indel(localAlign);
    ins = indels@insertion@unlistData;     
    del = indels@deletion@unlistData;
    
    snpsum = sum(snps$subject$Count);
    inssum = sum(ins@width)
    delsum = sum(del@width)
    
    tm = snpsum + inssum + delsum;
    
    nMuts = c(nMuts, tm);
  }
  
  predictedr = round(median(nMuts),0);
  predictedr = min(predictedr,opnuMutp)
  
  n = diversityRegionLength;
  
  if(opnuMut/diversityRegionLength > 0.9999)
  {
    n = 1.5*max(opnuMut,diversityRegionLength)
  }
  createPdfMat(n,predictedr,opnuMut)
  
  expMat = read.table("expV.txt")
  MMpredictedNop = mean(which(abs(expMat-opnuMut)==min(abs(expMat-opnuMut))))
  MMpredictedNopp = mean(which(abs(expMat-opnuMutp)==min(abs(expMat-opnuMutp))))
  
  system("rm expV.txt")
  
  write(x = c(paste("Ns = ",as.character(MMpredictedNop),sep=""),
              paste("Np = ",as.character(MMpredictedNopp),sep=""),
              paste("f_min = ",as.character(cutoff),sep="")),
        file = "ViQuaS-Richness.txt",ncolumns = 1,append = F,sep = "\t")
}

createPdfMat<-function(n,r,ipExpV)
{
  u = c(0:n);
  ncols = length(u);
  
  pdfMat = rep.int(x = 0, times = ncols)
  
  pdfMat[(r+1)] = 1;
  expV = c(sum(pdfMat*u));
  t=1;
  
  while(expV[length(expV)] < ipExpV)
  {
    t = t+1;
    
    pdfMatNew = rep.int(x = 0, times = ncols)
    
    for(v in u)
    {
      if((v>=r)&&(v<=(t*r)))
      {   
        p = 0;
        
        for(i in 0:r)
        {
          np = choose(n-v+i,i)*choose(v-i,r-i)/choose(n,r); 
          p = p + pdfMat[v-i+1]*np
        }
        pdfMatNew[(v+1)] = p;
      }
    }
    
    pdfMat = pdfMatNew;
    expV = c(expV, sum(pdfMat*u));
  }
  
  write.table(x=expV,file=paste("expV.txt",sep=""),append=F,row.names=F,col.names=F,sep="\t")
  
}
