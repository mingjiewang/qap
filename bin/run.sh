##filtration
#perl RawDataFiltration.pl -1 /home/wmj/projects/yzt/data/Y01_R1.fq.gz -2 /home/wmj/projects/yzt/data/Y01_Y01_R2.fil.fq.gz -l 250 -q 30 -o yzt/Y01

##barcode spliter
#perl BarcodeSplitter.pl -1 yzt/Y01/1.QualityFilter/Y01_R1.fil.fq -2 yzt/Y01/1.QualityFilter/Y01_Y01_R2.fil.fil.fq -b /home/wmj/projects/yzt/doc/barcode.txt -m 1 -r 1 -p BOS -o yzt/Y01/2.BarcodeSplit

##extract read id
#perl ExtractReadID.pl \
	#	-1	yzt/Y01/2.BarcodeSplit/P01_R1.fq,yzt/Y01/2.BarcodeSplit/P02_R1.fq,yzt/Y01/2.BarcodeSplit/P03_R1.fq,yzt/Y01/2.BarcodeSplit/P04_R1.fq,yzt/Y01/2.BarcodeSplit/P05_R1.fq,yzt/Y01/2.BarcodeSplit/P06_R1.fq,yzt/Y01/2.BarcodeSplit/P07_R1.fq,yzt/Y01/2.BarcodeSplit/P08_R1.fq,yzt/Y01/2.BarcodeSplit/P09_R1.fq,yzt/Y01/2.BarcodeSplit/P10_R1.fq \
	#	-2	yzt/Y01/2.BarcodeSplit/P01_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P02_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P03_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P04_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P05_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P06_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P07_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P08_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P09_Y01_R2.fil.fq,yzt/Y01/2.BarcodeSplit/P10_Y01_R2.fil.fq \
	#	-l	P01,P02,P03,P04,P05,P06,P07,P08,P09,P10 \
	#	-t  20 \
	#	-o  yzt/Y01/3.ReadID/
	
##extract seq
#for patient in P01 P02 P03 P04 P05 P06 P07 P08 P09 P10
#do
	#	perl ExtractSeqByIDInR.pl -1 yzt/Y01/1.QualityFilter/Y01_R1.fil.fq -2 yzt/Y01/1.QualityFilter/Y01_Y01_R2.fil.fil.fq -i yzt/Y01/3.ReadID/${patient}.readID.txt -f fastq -s $patient -o yzt/Y01/4.Sequence/
#done

##mapping
perl MapReadsToRef.pl \
		-1 yzt/Y01/4.Sequence/P01_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P02_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P03_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P04_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P05_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P06_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P07_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P08_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P09_Y01_R1.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P10_Y01_R1.fil.extractedSeq.fq \
		-2 yzt/Y01/4.Sequence/P01_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P02_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P03_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P04_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P05_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P06_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P07_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P08_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P09_Y01_R2.fil.extractedSeq.fq,yzt/Y01/4.Sequence/P10_Y01_R2.fil.extractedSeq.fq \
		-f sam \
		-t 20 \
		-p bowtie2 \
		-s pos \
		-l P01,P02,P03,P04,P05,P06,P07,P08,P09,P10 \
		-r ./test/ref/hbv.fasta \
		-o yzt/Y01/5.Map
	
	
		