import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class Base
{
	public double position;
	public char reference;
	public char consensus;
	public double A;
	public double C;
	public double G;
	public double T;
	public double del;
	public double coverage;
	public double probCoverage;
	public boolean corrected;
	public boolean removed;
	public double entropy;
	
	public Base()
	{
		position=-1;
		reference='X';
		consensus='X';
		A=0;
		C=0;
		G=0;
		T=0;
		del=0;
		coverage=0;
		probCoverage=-1;
		corrected=false;
		removed=false;
		entropy=-1;
	}
	
	public Base(double p, char r, char c, double b1, double b2, double b3, double b4, double b5, double cv, double ent)
	{
		position=p;
		reference=r;
		consensus=c;
		A=b1;
		C=b2;
		G=b3;
		T=b4;
		del=b5;
		coverage=cv;
		probCoverage=-1;
		corrected=false;
		removed=false;
		entropy=ent;
	}
	
	public void calculateEntropy()
	{
		if (coverage==0)
			entropy=-1;
		else
		{
			double ent=0;
			if (A>0)
				ent+=A/coverage*Math.log(A/coverage)/Math.log(2);
			if (C>0)
				ent+=C/coverage*Math.log(C/coverage)/Math.log(2);
			if (G>0)
				ent+=G/coverage*Math.log(G/coverage)/Math.log(2);
			if (T>0)
				ent+=T/coverage*Math.log(T/coverage)/Math.log(2);
			if (del>0)
				ent+=del/coverage*Math.log(del/coverage)/Math.log(2);
			if (ent<0)
				entropy=-ent;
			else
				entropy=0;
		}
	}
	
	public void setConsensus(double avg, double std, double threshold)
	{
		if (coverage<(A + C + G + T + del))
			coverage=(A + C + G + T + del);
		double refFreq = coverage - (A + C + G + T + del);
		if (reference=='A')
			A=refFreq;
		if (reference=='C')
			C=refFreq;
		if (reference=='G')
			G=refFreq;
		if (reference=='T')
			T=refFreq;
		if (reference=='-')
			del=refFreq;
		consensus='A';
		double maxFreq=A;
		if (C>maxFreq)
		{
			consensus='C';
			maxFreq=C;
		}
		if (G>maxFreq)
		{
			consensus='G';
			maxFreq=G;
		}
		if (T>maxFreq)
		{
			consensus='T';
			maxFreq=T;
		}
		if (del>maxFreq)
		{
			consensus='-';
			maxFreq=del;
		}
		probCoverage = Functions.zetaStandardProbability(((coverage-avg)/std));
		if (coverage<30 || probCoverage<threshold)
		{
			consensus='?';
			removed=true;
		}
	}
	
	public void correct(double tolerance, double erNoHomopol, double erHomopol, String genome, double mt)
	{
		if (consensus=='?')
			return;
		double errorRate = erNoHomopol;
		boolean indel = false;
		if (reference=='-' || consensus=='-' || Math.floor(position)!=position)
			indel = true;
		if (Functions.isHomopolymeric((int)Math.floor(position-1),genome,indel))
			errorRate = erHomopol;
		
		double pErrA = Functions.poissonError(A,coverage,errorRate);
		double pErrC = Functions.poissonError(C,coverage,errorRate);
		double pErrG = Functions.poissonError(G,coverage,errorRate);
		double pErrT = Functions.poissonError(T,coverage,errorRate);
		double pErrDel = Functions.poissonError(del,coverage,errorRate);
		
		pErrA = Math.min(1,pErrA*mt);
		pErrC = Math.min(1,pErrC*mt);
		pErrG = Math.min(1,pErrG*mt);
		pErrT = Math.min(1,pErrT*mt);
		pErrDel = Math.min(1,pErrDel*mt);
		
		boolean cazz=false;
		double remainder=0;
		if (A>0 && pErrA>tolerance)
		{
			corrected=true;
			remainder+=A;
			A=0;
		}
		if (C>0 && pErrC>tolerance)
		{
			corrected=true;
			remainder+=C;
			C=0;
		}
		if (G>0 && pErrG>tolerance)
		{
			corrected=true;
			remainder+=G;
			G=0;
		}
		if (T>0 && pErrT>tolerance)
		{
			corrected=true;
			remainder+=T;
			T=0;
		}
		if (del>0 && pErrDel>tolerance)
		{
			corrected=true;
			remainder+=del;
			del=0;
		}
		if (consensus=='A')
			A+=remainder;
		if (consensus=='C')
			C+=remainder;
		if (consensus=='G')
			G+=remainder;
		if (consensus=='T')
			T+=remainder;
		if (consensus=='-')
			del+=remainder;
		if (consensus=='-' && ((A+C+G+T)==0) && reference=='-')
			removed=true;
	}
}
