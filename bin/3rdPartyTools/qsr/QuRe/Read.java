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

public class Read
{
	public int idx;
	public String name;
	public String sequence;
	public String orientation;
	public String SNP_string;
	public Hashtable<String,String> SNP_hash;
	public Vector<SNP> SNP_list;
	public int mappingPosition;
	public float start;
	public float stop;
	public float frequency;
	public float coverage;
	public float prevalence;
	public float score;
	public float similarity;
	public int insertions;
	public double pvalue;
	public double adjustedPvalue;
	
	public Read()
	{
		idx=-1;
		name=null;
		sequence=null;
		SNP_string=null;
		SNP_hash=null;
		SNP_list=null;
		mappingPosition=-1;
		start=-1;
		stop=-1;
		frequency=-1;
		coverage=-1;
		prevalence=-1;
		score=-1;
		similarity=-1;
		insertions=-1;
		pvalue=-1;
		adjustedPvalue=-1;
	}
	
	public void correct(HashMap<Double,Base> baseSet)
	{
		Read r = this;
		Hashtable<String,String> newSNPHash = new Hashtable();
		Object [] key = baseSet.keySet().toArray();
		for (int i=0; i<key.length; i++)
		{
			Base b = baseSet.get(key[i]);
			if (r.SNP_hash.get(b.reference+"_"+b.position)!=null)
			{
				if (r.SNP_hash.get(b.reference+"_"+b.position).equals("A") && b.A>0)
					newSNPHash.put(b.reference+"_"+b.position,"A");
				if (r.SNP_hash.get(b.reference+"_"+b.position).equals("C") && b.C>0)
					newSNPHash.put(b.reference+"_"+b.position,"C");
				if (r.SNP_hash.get(b.reference+"_"+b.position).equals("G") && b.G>0)
					newSNPHash.put(b.reference+"_"+b.position,"G");
				if (r.SNP_hash.get(b.reference+"_"+b.position).equals("T") && b.T>0)
					newSNPHash.put(b.reference+"_"+b.position,"T");
				if (r.SNP_hash.get(b.reference+"_"+b.position).equals("-") && b.del>0)
					newSNPHash.put(b.reference+"_"+b.position,"-");
			}
			if 
			(
				r.SNP_hash.get(b.reference+"_"+b.position)==null && 
				(r.start<=b.position) && (r.stop>=b.position) &&
				(b.reference!=b.consensus) &&
				( (b.A==0 && b.reference=='A') || (b.C==0 && b.reference=='C') || (b.G==0 && b.reference=='G') || (b.T==0 && b.reference=='T') || (b.del==0 && b.reference=='-') )
			)
				newSNPHash.put(b.reference+"_"+b.position,b.consensus+"");
		}
		r.SNP_hash=newSNPHash;
		r.updateSNPFromHash();
	}
	
	public void setPvalue(double rndScoreAvg, double rndScoreStd)
	{
		double zeta = (this.score-rndScoreAvg)/rndScoreStd;
		this.pvalue = Functions.alignmentScoreSignificance(zeta);
	}
	
	public void align(Hashtable genomeDictionary, String refGenome, double avgReadLength, double stdReadLength, Matrix matrix, float gop, float gep, int kmer)
	{
		Logger l = Logger.getLogger(SmithWatermanGotoh.class.getName());
		l.setLevel(Level.OFF);
		int bound = (int)(avgReadLength+3*stdReadLength);
		Integer [] fwdSBLAT = Functions.SBLAT(genomeDictionary,this.sequence,kmer);
		int fwdMapPos = Functions.getMappingPosition(fwdSBLAT);
		int fwdStart = Math.max(fwdMapPos-bound,0);
		int fwdStop = Math.min(fwdMapPos+bound,refGenome.length());
		String fwdGenomeCut = refGenome.substring(fwdStart,fwdStop);
		if (fwdMapPos==-1)
			fwdGenomeCut = refGenome;					
		Sequence fwdQuery = new Sequence(this.sequence,"","",0);
		Sequence fwdReference = new Sequence(fwdGenomeCut,"","",0);
		Alignment fwdAlignment = new Alignment();
		fwdAlignment = SmithWatermanGotoh.align(fwdQuery, fwdReference, matrix, gop, gep);
		float fwdScore = fwdAlignment.getScore();
			
		Integer [] rwdSBLAT = Functions.SBLAT(genomeDictionary,Functions.reverseComplement(this.sequence),kmer);
		int rwdMapPos = Functions.getMappingPosition(rwdSBLAT);
		int rwdStart = Math.max(rwdMapPos-bound,0);
		int rwdStop = Math.min(rwdMapPos+bound,refGenome.length());
		String rwdGenomeCut = refGenome.substring(rwdStart,rwdStop);
		if (rwdMapPos==-1)
			rwdGenomeCut = refGenome;					
		Sequence rwdQuery = new Sequence(Functions.reverseComplement(this.sequence),"","",0);
		Sequence rwdReference = new Sequence(rwdGenomeCut,"","",0);
		Alignment rwdAlignment = new Alignment();
		rwdAlignment = SmithWatermanGotoh.align(rwdQuery, rwdReference, matrix, gop, gep);
		float rwdScore = rwdAlignment.getScore();

		Alignment alignment = fwdAlignment;
		this.orientation="forward";
		this.mappingPosition=fwdMapPos;
		if (rwdScore>fwdScore)
		{
			alignment = rwdAlignment;
			this.orientation="reverse";
			this.mappingPosition=rwdMapPos;
		}
		this.similarity=(float)(alignment.getSimilarity())/(float)(alignment.getSequence1().length);
		this.start=1+alignment.getStart2()+Math.max(this.mappingPosition-bound,0);
		this.stop=this.start+alignment.getSequence1().length;
		this.score=alignment.getScore();
		this.setSNP(alignment);
		this.sequence=String.valueOf(alignment.getSequence1()).replaceAll("-","");
	}
	
	public void setSNP(Alignment a)
	{
		Vector<SNP> sl = new Vector();
		Hashtable sh = new Hashtable();
		String ss = "";
		char [] query = a.getSequence1();
		char [] refer = a.getSequence2();
		double ins = 0;
		//double frac = 0.5f;
		double frac = 0.000000001f;
		double insCount = 0;
		for (int i=0; i<query.length; i++)
		{
			double pos = (this.start+i)-insCount;
			if (query[i]!=refer[i] && refer[i]!='-')
			{
				SNP snp = new SNP(pos,refer[i],query[i]);
				sl.add(snp);
				sh.put(refer[i]+"_"+pos,query[i]+"");
				ss=ss+refer[i]+"_"+pos+"_"+query[i]+",";
				ins=pos;
				//frac=0.5f;
				frac=0.000000001f;
			}
			else
				if (query[i]!=refer[i] && refer[i]=='-')
				{
					ins=ins+frac;
					SNP snp = new SNP(ins,refer[i],query[i]);
					sl.add(snp);
					sh.put(refer[i]+"_"+ins,query[i]+"");
					ss=ss+refer[i]+"_"+ins+"_"+query[i]+",";
					insCount++;
					//frac=frac/2f;
					frac=frac+0.000000001f;
				}
				else
				{
					ins=pos;
					//frac=0.5f;
					frac=0.000000001f;
				}
		}
		SNP_list=sl;
		SNP_hash=sh;
		SNP_string=ss;
		stop=start+a.getSequence1().length-(float)(insCount);
		insertions=(int)(insCount);
	}
	
	public void updateSNPFromHash()
	{
		String newSNPString = "";
		Vector<SNP> newSNPList = new Vector();
		Object [] key = this.SNP_hash.keySet().toArray();
		for   (int i=0; i<key.length; i++)
		{
			String s = this.SNP_hash.get(key[i]);
			char ref = key[i].toString().split("_")[0].charAt(0);
			double pos = Double.parseDouble(key[i].toString().split("_")[1]);
			char mut = s.charAt(0);
			SNP newSNP = new SNP(pos,ref,mut);
			newSNPList.add(newSNP);
		}
		this.SNP_list = newSNPList;
		Comparator<SNP> comparatorPositions = new Comparator<SNP>()
		{
			public int compare(SNP r1, SNP r2)
			{
				if (r1.position<r2.position)
					return -1;
				if (r1.position>r2.position)
					return 1;
				return 0;
			}
		};
		Collections.sort(this.SNP_list,comparatorPositions);
		for (int i=0; i<this.SNP_list.size(); i++)
		{
			SNP snp = this.SNP_list.get(i);
			newSNPString = newSNPString + snp.consensus + "_" + snp.position + "_" + snp.base + ",";
		}
		this.SNP_string = newSNPString;
	}
	
	public String getSNPString(double sta, double sto)
	{
		String thisSNP="";
		for (int i=0; i<this.SNP_list.size(); i++)
		{
			SNP snp = this.SNP_list.get(i);
			if (snp.position>=sta && snp.position<=sto)
				thisSNP+=snp.toString();
		}
		return thisSNP;
	}
	
	public double distance(Read r, double sta, double sto)
	{
		double distance = 0;
		String s1 = this.getSNPString(sta,sto);
		String s2 = r.getSNPString(sta,sto);
		String [] SNP1 = s1.split(",");
		String [] SNP2 = s2.split(",");
		LinkedList l1 = new LinkedList(Arrays.asList(SNP1));
		LinkedList l2 = new LinkedList(Arrays.asList(SNP2));
		for (int i=0; i<SNP1.length; i++)
			if (!l2.contains(SNP1[i]))
				distance++;
		for (int i=0; i<SNP2.length; i++)
			if (!l1.contains(SNP2[i]))
				distance++;
		//System.out.println(s1);
		//System.out.println(s2);
		//System.out.println(distance);
		return distance;
	}
	
	public boolean spans(double sta, double sto)
	{
		if (start<=sta && stop>=sto)
			return true;
		return false;
	}
	
	public boolean overlaps(Read r)
	{
		if ((start<r.start && stop<r.stop && stop>r.start) || (r.start<start && r.stop<stop && r.stop>start))
			return true;
		return false;
	}

	public double locationDistance(Read r)
	{
		return Math.sqrt(Math.pow((this.start-r.start),2)+Math.pow((this.stop-r.stop),2));
	}
}
