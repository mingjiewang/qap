import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class LocalVariant
{
	public String amplSNP;
	public String overlSNP1;
	public String overlSNP2;
	public String SNP;
	public double frequency;
	
	public LocalVariant()
	{
		amplSNP=null;
		overlSNP1=null;
		overlSNP2=null;
		SNP=null;
		frequency=-1;
	}
	public LocalVariant(String asnp, String osnp1, String osnp2, String snp, double f)
	{
		amplSNP=asnp;
		overlSNP1=osnp1;
		overlSNP2=osnp2;
		SNP=snp;
		frequency=f;
	}
	public boolean overlaps(LocalVariant b)
	{
		return (this.overlSNP2.equals(b.overlSNP1));
	}
	public int searchMateFwd(LocalVariantSet lvs)
	{
		for(int j=0; j<lvs.lvsA.size(); j++)
			if (this.overlaps(lvs.lvsA.get(j)))
				return j;
		return -1;
	}
	public int searchMateRwd(LocalVariantSet lvs)
	{
		for(int j=0; j<lvs.lvsA.size(); j++)
			if (lvs.lvsA.get(j).overlaps(this))
				return j;
		return -1;
	}
	public double distance(LocalVariant b)
	{
		double distance = 0;
		String [] SNP1 = this.SNP.split(",");
		String [] SNP2 = b.SNP.split(",");
		LinkedList l1 = new LinkedList(Arrays.asList(SNP1));
		LinkedList l2 = new LinkedList(Arrays.asList(SNP2));
		for (int i=0; i<SNP1.length; i++)
			if (!l2.contains(SNP1[i]))
				distance++;
		for (int i=0; i<SNP2.length; i++)
			if (!l1.contains(SNP2[i]))
				distance++;
		return distance;
	}
	public String snpDifferences(LocalVariant b)
	{
		String diff = "";
		String [] SNP1 = this.SNP.split(",");
		String [] SNP2 = b.SNP.split(",");
		LinkedList l1 = new LinkedList(Arrays.asList(SNP1));
		LinkedList l2 = new LinkedList(Arrays.asList(SNP2));
		for (int i=0; i<SNP1.length; i++)
			if (!l2.contains(SNP1[i]))
				diff+=SNP1[i]+",";
		for (int i=0; i<SNP2.length; i++)
			if (!l1.contains(SNP2[i]))
				diff+=SNP2[i]+",";
		return diff;
	}
	public double lengthUpdate(LocalVariant b)
	{
		String [] SNP1 = this.SNP.split(",");
		String [] SNP2 = b.SNP.split(",");
		LinkedList l1 = new LinkedList(Arrays.asList(SNP1));
		LinkedList l2 = new LinkedList(Arrays.asList(SNP2));
		LinkedList<String> dels = new LinkedList();
		LinkedList<String> inss = new LinkedList();
		
		for (int i=0; i<SNP1.length; i++)
			if (SNP1[i].indexOf("_-")!=-1)
				if (!dels.contains(SNP1[i]))
					dels.add(SNP1[i]);
		for (int i=0; i<SNP2.length; i++)
			if (SNP2[i].indexOf("_-")!=-1)
				if (!dels.contains(SNP2[i]))
					dels.add(SNP2[i]);
		for (int i=0; i<SNP1.length; i++)
			if (SNP1[i].indexOf("-_")!=-1)
				if (!inss.contains(SNP1[i]))
					inss.add(SNP1[i]);
		for (int i=0; i<SNP2.length; i++)
			if (SNP2[i].indexOf("-_")!=-1)
				if (!inss.contains(SNP2[i]))
					inss.add(SNP2[i]);
		double plus = inss.size();
		double minus = 0;
		for (int i=0; i<dels.size(); i++)
		{
			String s = dels.get(i);
			if (this.SNP.indexOf(s)!=-1 && b.SNP.indexOf(s)!=-1)
				minus++;
		}
		double res = plus - minus;
		return res;
	}
	public int classify(Vector<LocalVariant> clusterSet)
	{
		int ind = -1;
		double dist = Double.MAX_VALUE;
		for (int i=0; i<clusterSet.size(); i++)
		{
			double d = this.distance(clusterSet.get(i));
			if (d<dist)
			{
				ind=i;
				dist=d;
			}
		}
		return ind;
	}
}
