LICENSE
PredictHaplo-Paired is free software.  The text file COPYING, located in this directory, contains copyright and licensing details. The software package contains a copy of the SCYTHE C++ library (http://scythe.wustl.edu/) which is also published under the GPL.

GENERAL
This software aims at reconstructing haplotypes from next-generation sequencing data. A typical application example is to reconstruct HIV-1 haplotypes present in a blood sample from a patient based on paired-end Illumina reads. 


INSTALLATION
PredictHaplo has been tested with gcc 4.6.3 under Ubuntu Linux. 
It should compile on similar platforms with gcc version 4.4.7 or greater.

1) Create a directory with absolute path <yourPATH>

       mkdir <yourPATH>

2) Copy "PredictHaplo-Paired-0.4.tgz" to this new directory.
3) Inside this directory, unpack the "PredictHaplo-Paired-0.4.tgz" tar-bundle:

       tar xfvz PredictHaplo-Paired-0.4.tgz
       
4) Change to the "<yourPATH>/PredictHaplo-Paired-0.4"  directory. Unpack the SCYTHE library tar-bundle:
       
       tar xfvz scythestat-1.0.3.tar.gz


5) Change to the "<yourPATH>/PredictHaplo-1.0/scythestat-1.0.3" directory and install SCYTHE:

       ./configure --prefix=<yourPATH>/PredictHaplo-Paired-0.4/NEWSCYTHE
       make install

6) Go back to "<yourPATH>/PredictHaplo-Paired-0.4/" and install PredictHaplo-Paired:

      make

7) Done, hopefully :-)


DOCUMENTATION
Here are some basic informations:

The program is invoked with a configuration file as "PredictHaplo-Paired <config_file>", see the hopefully self-contained example file "config_5V". The syntax of this config file is quite simple: lines starting with "%" are treated as comments. 

The first non-comment line specifies a prefix for all files that are created during execution. 
The second such line contains the filename of the reference sequence (in FASTA format), the third line specifies the way how intermediate (local) reconstructions are visualized for better inspection (2: in every local window a html file with all reconstruction results and all aligned reads in this window is stored, 1: only the "interesting" positions are shown in the html file, where interesting mans that they pass an entropy threshold, 0: no visualization).

Then comes the filename of the aligned reads (SAM format), the reads must be aligned to the reference sequence provided in the first line. One way to create such a SAM file is to use the "bwa" aligner (http://bio-bwa.sourceforge.net/) with commands like:
    bwa index <Reference.fasta>
    bwa mem <Reference.fasta> <Reads.fastq>  >  alignedReads.sam
Then comes a flag indicating if there are known "true" haplotypes, followed by the name of the corresponding FASTA file that contains the aligned true haplotypes (alignment w.r.t. the same reference sequence as above. The alignment must have the same length as the reference sequence, ignore all positions with insertions into the reference!). Fill in some arbitrary string if you don't have such ground-truth. 

The next line is a flag indicating if a local analysis should be performed (must be 1 when run the first time), then comes the maximum number of reads selected in one of the local windows (3000-20000 should be convenient), and finally some entropy threshold which basically specifies the smallest expected haplotype frequency. 

Several additional parameters can be specified: start + stop position of the reconstruction (with respect to the provided reference genome), a quality threshold parameter that selects reads with good alignments to the reference genome, the minimum read length, the maximum fraction of gaps within a read (relative to alignment length), the minimum ratio of alignment score to alignment length, the prior parameter of the multinomial probability tables of the nucleotides at every position, the "min_overlap_factor" (reads must have an overlap with the local reconstruction window of at least this factor times the window size), the "local_window_size_factor" (size of  local reconstruction window relative to average read length), the maximum number of clusters in the truncated Dirichlet process, and the number of MCMC iterations. The last line is an indicator variable for handling deletions in the reads (0 means ignoring deletions; this is useful in regions where one does not expect "true" deletions).

DEMO
The package includes the Illumina reads from the "5-virus-mix" dataset that can be found in the supplement of [Di Giallonardo, Töpfer et al. Full-length haplotype reconstruction to infer the structure of heterogeneous virus populations. NAR, 2014 10.1093/nar/gku537], aligned to the hxb2 reference "5V_ref_seq.fasta", and the "true" haplotypes aligned to the same reference: "5VirusMixReference.fasta".  With the config file provided, you should be able to successfully detect the 5 haplotypes in the gag-pol region.  During execution, progress can be monitored by looking at the created html-files (one per local or "global" window). These html files contain the reference sequence (move with the mouse pointer over a letter to see the position). If true haplotypes are available, these are shown as well, and also the "true" mutations are highlighted in the reference sequence. Then follow the reconstructed haplotypes (in the current analysis window), together with some information about their quality: the "Overlaps" section, where an entry "10:20" means that in every position in the window there are at least 10 reads assigned to one haplotype which have an overlap with this position of 20 (only positions that pass the entropy criterion are counted), and that there is one (or more) position where there are only 10 such reads. Any number x in "10:x" smaller then 15 might indicate that there are problems with this haplotype. If true haplotypes are available, the best matching ones are shown, and the number of errors ("cost"), and the position of these errors (if there are any).
The files for the "local" windows additionally contain all the aligned reads ordered by haplotype assignment. The inferred haplotypes are also output in fasta format (filenames have the following structure: prefix + local/global + window_start + window_stop + .fas).

