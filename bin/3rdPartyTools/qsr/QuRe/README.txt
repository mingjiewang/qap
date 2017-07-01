QuRe (v0.99971) is distributed under the GNU general public license. It has been released only for research, academic and non-profit purposes. The software should not be used as a clinical decision tool.

COMMAND LINE USAGE SYNTAX:

"java [-classpath .] [-Xmx{1,2,3,...}G] QuRe read_file reference_genome_file [homopolymericErrorRate nonHomopolymericErrorRate iterations]"

The read file and the reference genome file must be in FASTA format (nucleotides). The read file can contain ambiguous base codes (i.e. R/Y/K/M/S/W/B/D/H/V/N). The reference file DOES NOT ALLOW for ambiguous base codes, and all non-ACGT characters are ignored.

If the last three parameters are not inserted, default values are used (0.01, 0.005, 3000).

The user might get a "NoclassDefFoundError" if the classpath is not setup or referenced correctly. In case, run the program with the -classpath . option ("." is the local directory), referencing (if needed) also the path and sub-paths to java, QuRe and the internal QuRe subdirectories.

Note that the provided class files have been compiled using using the 64 bit Java Development Kit, Standard Edition 6 (http://www.oracle.com/technetwork/java), ), coupled with the pre-compiled open-source Java library JAligner (http://jaligner.sourceforge.net/). If the user wants to compile again the QuRe code, use simply "javac *.java" from the command line.

The option "-Xmx..." can be used to increase the default memory usage of the Java Virtual Machine.

MEMORY REQUIREMENTS:

For processing a fasta file of ~16,500 reads (avg. read length of ~350 bases) with a reference genome of ~10,000 bases, and standard parameter settings, QuRe requires ~1Gb of RAM. For a file of ~20,000 or ~45,000 reads and similar configuration/characteristics as above, the program usually assesses to ~1.5Gb of maximal RAM usage.

RUNNING WITHOUT THE ERROR CORRECTION MODULE

QuRe can be run overriding the error correction module (e.g. if reads were corrected with another program) by setting up the error parameters very close to zero (say 1E-25), for instance "java -Xmx3G QuRe read_file reference_genome_file 1E-25 1E-25 1000". Not sure if numerical errors may arise with exact zero values.

OUTPUT FILES:

"filename_alignedReads.txt"
	This file reports reads that align significantly to the reference genome, printing the aligned sequence (trimmed and gap-stripped), the corresponding start/stop positions of the pairwise alignment and the p-value of the alignment. Base changes from the reference are reported in the format r_p_b, i.e. reference base (r), reference position (p), and replaced base (b). Insertions and deletions are encoded with "-". Insertions and multiple insertions take a fractional position. For instance a C insertion at position 150 is encoded as -_150.5_C, and any other subsequent insertion is added a half of the previous step (i.e. 150.75, 150.875, ...). Base changes are corrected according to the Poisson-based error model.
"filename_snpTable.txt"
	This file prints the positions with respect to the reference genome numbering that were significantly covered by the read set, with the corresponding consensus base and the specific base composition for each position. Read coverage and base composition entropy are also reported. The file can be used to look at the genome variation base by base singularly.
"filename_overlappingIntervalsSet.txt"
	This file reports the optimal set of overlapping intervals, with start, overlap and stop positions with respect to the reference genome. For each overlapping interval, all the distinct reads found are reported along with their relative frequency (%). Reads are defined by their list of base changes from the reference genome, divided into changes in the overlapping and non-overlapping parts. Each overlapping interval can be regarded as a local quasispecies reconstruction, i.e. the distinct variants found in a particular sub-region of the mapped reference genome.
"filename_reconstructedVariantsNoClustering.txt"
	This is an INTERMEDIATE fasta file where reconstructed variants (i.e. global quasispecies reconstruction across the whole, mapped reference genome) are reported along with their estimated prevalence (%) BEFORE THE FINAL CLUSTERING IS PERFORMED.
"filename_reconstructedVariants.txt"
	This is the final fasta file where reconstructed variants (i.e. global quasispecies reconstruction across the whole, mapped reference genome) are reported along with their estimated prevalence (%).

Please refer to the paper cited below for a detailed description of methods.

CITATION:

Please cite the following when using QuRe in publications or research: "Prosperi MC, and Salemi, M. QuRe: software for viral quasispecies characterization from next-generation sequencing data. Bioinformatics. 2012 Jan 1;28(1):132-3".


------------------------------------------------------------
                    CONTACT INFORMATION
------------------------------------------------------------
Mattia Prosperi, Ph.D.

University of Manchester, UK

e-mail: ahnven@gmail.com

QuRe is no longer supported by grants, so (small) donations
via PayPal are kindly accepted: ahnven@gmail.com
------------------------------------------------------------