import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class GlobalVariant
{
	public String SNP;
	public String sequence;
	public double frequency;
	public double stdevFreq;
	
	public GlobalVariant()
	{
		SNP=null;
		sequence=null;
		frequency=-1;
		stdevFreq=-1;
	}
	public GlobalVariant(String snp, String seq, double f, double s)
	{
		SNP=snp;
		sequence=seq;
		frequency=f;
		stdevFreq=s;
	}
	public void setSequence(String refGenome, double start, double stop)
	{
		HashMap<Double,Base> bs = new HashMap(refGenome.length());
		for (int i=(int)(start-1); i<(int)(stop); i++)
		{
			Base b =  new Base();
			b.reference=refGenome.charAt(i);
			b.position=(double)(i+1);
			bs.put(b.position,b);
		}
		if (SNP!=null && SNP.length()>0)
		{
			if (SNP.startsWith(","))
				SNP=SNP.substring(1);
			if (!SNP.equals(""))
			{
				String [] snpA = SNP.split(",");
				for (int i=0; i<snpA.length; i++)
				{
					String [] snp = snpA[i].split("_");
					Base b =  new Base();
					b.reference=snp[2].charAt(0);
					b.position=Double.parseDouble(snp[1]);
					bs.put(b.position,b);
				}
			}
		}
		String s="";
		Object [] key = bs.keySet().toArray();
		Arrays.sort(key);
		for   (int i=0; i<key.length; i++)
		{
			char c = bs.get(key[i]).reference;
			s+=c;
		}
		sequence=s.replaceAll("-","");
	}
}
