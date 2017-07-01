import java.io.*;
import java.lang.*;
import java.util.*;
import java.util.logging.*;
import java.math.*;
import jaligner.Alignment;
import jaligner.Sequence;
import jaligner.SmithWatermanGotoh;
import jaligner.formats.Pair;
import jaligner.formats.CLUSTAL;
import jaligner.formats.FASTA;
import jaligner.matrix.*;
import jaligner.util.SequenceParser;

public class SNP
{
	public double position;
	public char consensus;
	public char base;
	
	public SNP()
	{
		position=-1;
		consensus='X';
		base='X';
	}
	
	public SNP(double p, char c, char b)
	{
		position=p;
		consensus=c;
		base=b;
	}
	
	public String toString()
	{
		String s = consensus+"_"+position+"_"+base+",";
		return s;
	}
	
	public boolean equals(SNP b)
	{
		return (this.position==b.position && this.consensus==b.consensus && this.base==b.base);
	}
}
