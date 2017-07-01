import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class LocalVariantSetEnsemble
{
	public LocalVariantSet [] lvse;
	public double [] starts;
	public double [] stops;
	public double [] overl1;
	public double [] overl2;
	
	public LocalVariantSetEnsemble()
	{
		lvse=null;
		starts=null;
		stops=null;
		overl1=null;
		overl2=null;
	}
	
	public LocalVariantSetEnsemble(Vector<Read> p, AmpliconSet a, String genome, double hErr, double nonhErr, int n_proc)
	{
		starts=new double[a.starts.length];
		stops=new double[a.starts.length];
		overl1=new double[a.starts.length];
		overl2=new double[a.starts.length];
		lvse=new LocalVariantSet[a.starts.length];
		lvse[0]=new LocalVariantSet(p,a.starts[0],a.starts[0],a.starts[1],a.stops[0],genome,hErr,nonhErr,n_proc);
		starts[0]=a.starts[0];
		stops[0]=a.stops[0];
		overl1[0]=a.starts[0];
		overl2[0]=a.starts[1];
		for (int i=1; i<lvse.length-1; i++)
		{
			double start=a.starts[i];
			double middle1=a.stops[i-1];
			double middle2=a.starts[i+1];
			double stop=a.stops[i];
			LocalVariantSet lvs = new LocalVariantSet(p,start,middle1,middle2,stop,genome,hErr,nonhErr,n_proc);
			lvse[i]=lvs;
			starts[i]=a.starts[i];
			stops[i]=a.stops[i];
			overl1[i]=middle1;
			overl2[i]=middle2;
		}
		lvse[lvse.length-1]=new LocalVariantSet(p,a.starts[lvse.length-1],a.stops[lvse.length-2],a.stops[lvse.length-1],a.stops[lvse.length-1],genome,hErr,nonhErr,n_proc);
		starts[lvse.length-1]=a.starts[lvse.length-1];
		stops[lvse.length-1]=a.stops[lvse.length-1];
		overl1[lvse.length-1]=a.stops[lvse.length-2];
		overl2[lvse.length-1]=a.stops[lvse.length-1];
	}
	
	public int getGuideDistribution()
	{
		int maxNumVariants=0;
		for (int i=0; i<lvse.length; i++)
		{
			LocalVariantSet lvs = lvse[i];
			int numVariants = lvs.lvsA.size();
			if (numVariants>maxNumVariants)
				maxNumVariants=numVariants;
		}
		double [][] multinomialDistributions = new double [lvse.length][maxNumVariants];
		for (int i=0; i<lvse.length; i++)
			for (int j=0; j<maxNumVariants; j++)
			{
				if (j<lvse[i].lvsA.size())
					multinomialDistributions[i][j]=lvse[i].lvsA.get(j).frequency;
				else
					multinomialDistributions[i][j]=0;
			}
		
		int index = (int)(Math.random()*(multinomialDistributions.length-1));
		double bestChiProb = -1*Double.MAX_VALUE;
		for (int i=0; i<multinomialDistributions.length; i++)
		{
			double chiProb = 0; 
			for (int j=0; j<multinomialDistributions.length; j++)
			{
				double chiSquare = 0;
				for (int k=0; k<lvse[i].lvsA.size(); k++)
					chiSquare+=(multinomialDistributions[i][k]-multinomialDistributions[j][k])*(multinomialDistributions[i][k]-multinomialDistributions[j][k])/multinomialDistributions[i][k];
				chiProb += Math.log(Functions.inverseChi(chiSquare,lvse[i].lvsA.size()-1));
			}
			if (chiProb>bestChiProb)
			{
				bestChiProb=chiProb;
				index=i;
			}
		}
		
		/*
		int index = (int)(Math.random()*(multinomialDistributions.length-1));
		double minSQM = Double.MAX_VALUE;
		for (int i=0; i<multinomialDistributions.length; i++)
		{
			double sqm = 0; 
			for (int j=0; j<multinomialDistributions.length; j++)
				for (int k=0; k<multinomialDistributions[0].length; k++)
					sqm+=(multinomialDistributions[i][k]-multinomialDistributions[j][k])*(multinomialDistributions[i][k]-multinomialDistributions[j][k]);
			if (sqm<minSQM)
			{
				minSQM=sqm;
				index=i;
			}
		}
		*/
		return index;
	}
	
	public void printToFile(String fileName)
	throws Exception
	{
		FileWriter fw = new FileWriter(fileName);
		fw.write("amplicon\tstart\toverl1\toverl2\tstop\toverlSNP1\tamplSNP\toverlSNP2\tSNP\tfrequency\r\n");
		for (int i=0; i<lvse.length; i++)
		{
			LocalVariantSet lvs = lvse[i];
			for (int j=0; j<lvs.lvsA.size(); j++)
			{
				LocalVariant lv= lvs.lvsA.get(j);
				fw.write(i+"\t");
				fw.write(Math.round(starts[i])+"\t");
				fw.write(Math.round(overl1[i])+"\t");
				fw.write(Math.round(overl2[i])+"\t");
				fw.write(Math.round(stops[i])+"\t");
				fw.write(lv.overlSNP1+"\t");
				fw.write(lv.amplSNP+"\t");
				fw.write(lv.overlSNP2+"\t");
				fw.write(lv.SNP+"\t");
				fw.write(lv.frequency+"\r\n");
			}	
		}
		fw.close();
	}
	
	public Vector<GlobalVariant> quasispeciesReconstructor(String refGenome)
	{
		//System.out.print("\texecuting core reconstruction algorithm ");
		
		System.out.print("\r\n\t\t");
			
		Vector<GlobalVariant> reconstructedQuasispecies = new Vector();
		while(true)
		{
			for (int i=0; i<(lvse.length-1); i++)
			{
				for (int j=0; j<lvse[i].lvsA.size(); j++)
				{
					boolean atLeastOneMatch = false;
					for (int k=0; k<lvse[i+1].lvsA.size(); k++)
					{
						if ( lvse[i].lvsA.get(j).overlaps(lvse[i+1].lvsA.get(k)) )
						{
							atLeastOneMatch = true;
							break;
						}
					}
					if (!atLeastOneMatch)
						lvse[i].lvsA.get(j).frequency=0;
				}
			
			}
		
			for (int i=0; i<lvse.length; i++)
			{
				//Vector<LocalVariant> lvA = lvse[i].lvsA;
				int j=0;
				while(j<lvse[i].lvsA.size())
				{
					if (lvse[i].lvsA.get(j).frequency<=0)
						lvse[i].lvsA.remove(j);
					else
						j++;
				}
			}
			
			double multinomialSize=this.lvse[0].lvsA.size();
			for (int k=1; k<lvse.length; k++)
			{
				if (this.lvse[k].lvsA.size()<multinomialSize)
					multinomialSize=this.lvse[k].lvsA.size();
			}
			if (multinomialSize==0)
				break;
			
			for (int i=0; i<lvse.length; i++)
				lvse[i].sortMultinomialFrequencies();
			
			int g=this.getGuideDistribution();
			//System.out.println("guide distribution is "+g);
			GlobalVariant gv = this.globalVariantReconstructor(g);
			
			double totComb=this.lvse[0].lvsA.size();
			multinomialSize=this.lvse[0].lvsA.size();
			for (int k=1; k<lvse.length; k++)
			{
				totComb=totComb*this.lvse[k].lvsA.size();
				if (this.lvse[k].lvsA.size()<multinomialSize)
					multinomialSize=this.lvse[k].lvsA.size();
			}
			
			//System.out.println("gv.frequency "+gv.frequency);
			if (gv.frequency==-1 || multinomialSize==0)
				break;
			else
			{
				//String perc = totComb+" putative comparisons left ";
				//System.out.print(perc);
				gv.setSequence(refGenome,starts[0],stops[stops.length-1]);
				reconstructedQuasispecies.add(gv);
				//for (int o=0; o<perc.length(); o++)
				//	System.out.print("\b");
			}
		}
		System.out.print("\b reconstruction(s) done.");
		//System.out.println();
		return reconstructedQuasispecies;
	}
	
	public void reconstruct(int i, int k, Vector<Integer> r, int g)
	{try
	{
		//System.out.print("i;k "+i+";"+k+"   dim="+r.size()+"   ");
		//for (int o=0; o<r.size(); o++)
		//	System.out.print(r.get(o)+" ");
		//System.out.println();
		if (i<g)
		{
			return;
		}
		if ( i>=(lvse.length-1) )
		{
			//System.out.println("!!!trovato!!!");
			return;
		}
		if ( k>=lvse[i+1].lvsA.size() )
		{
			int newK=r.get(i-g)+1;
			r.remove(i-g);
			this.reconstruct(i-1,newK,r,g);
			return;
		}
		if ( lvse[i].lvsA.get(r.get(i-g)).overlaps(lvse[i+1].lvsA.get(k)) )
		{
			r.add(k);
			this.reconstruct(i+1,0,r,g);
			return;
		}
		else
		{
			this.reconstruct(i,k+1,r,g);
			return;
		}
	}
	catch(Exception e)
	{
		//System.out.println("eccezione. "+e);
		System.exit(0);
	}}
	
	public void reconstructBack(int i, int k, Vector<Integer> r, int g)
	{try
	{
		//System.out.print("rwd i;k "+i+";"+k+"   dim="+r.size()+"   ");
		//for (int o=0; o<r.size(); o++)
		//	System.out.print(r.get(o)+" ");
		//System.out.println();
		if (i>g)
		{
			return;
		}
		if ( i<=0 )
		{
			//System.out.println("!!!trovatoRWD!!!");
			return;
		}
		if ( k>=lvse[i-1].lvsA.size() )
		{
			int newK=r.get(g-i)+1;
			r.remove(g-i);
			this.reconstructBack(i+1,newK,r,g);
			return;
		}
		if ( lvse[i-1].lvsA.get(k).overlaps(lvse[i].lvsA.get(r.get(g-i))) )
		{
			r.add(k);
			this.reconstructBack(i-1,0,r,g);
			return;
		}
		else
		{
			this.reconstructBack(i,k+1,r,g);
			return;
		}
	}
	catch(Exception e)
	{
		//System.out.println("eccezioneRwd. "+e);
		System.exit(0);
	}}
	
	public int getMismatchIndexFwd(int g, int [] indices)
	{
		for (int i=g; i<lvse.length-1; i++)
			if ( !lvse[i].lvsA.get(indices[i]).overlaps(lvse[i+1].lvsA.get(indices[i+1])) )
				return i+1;
		return -1;
	}
	
	public int getMismatchIndexRwd(int g, int [] indices)
	{
		for (int i=g; i>0; i--)
			if ( !lvse[i-1].lvsA.get(indices[i-1]).overlaps(lvse[i].lvsA.get(indices[i])) )
				return i-1;
		return -1;
	}
	
	public int [] updateMatchingSet(int g, int [] indices)
	{
		int [] updatedIndices = new int [lvse.length];
		for (int i=0; i<lvse.length; i++)
			updatedIndices[i]=indices[i];
		int mismatchIndexFwd = getMismatchIndexFwd(g, indices);
		int mismatchIndexRwd = getMismatchIndexRwd(g, indices);
		if (mismatchIndexFwd!=-1)
		{
			updatedIndices[mismatchIndexFwd]++;
			while (updatedIndices[mismatchIndexFwd]>=lvse[mismatchIndexFwd].lvsA.size() && mismatchIndexFwd>g)
			{
				mismatchIndexFwd--;
				updatedIndices[mismatchIndexFwd]++;
			}
			if (mismatchIndexFwd==g && updatedIndices[mismatchIndexFwd]>=lvse[mismatchIndexFwd].lvsA.size())
				return null;
			for (int i=mismatchIndexFwd+1; i<lvse.length; i++)
				updatedIndices[i]=0;
			if (mismatchIndexFwd==g)
			{
				for (int i=0; i<mismatchIndexFwd; i++)
					updatedIndices[i]=0;
				return updatedIndices;
			}
		}
		if (mismatchIndexRwd!=-1)
		{
			updatedIndices[mismatchIndexRwd]++;
			while (updatedIndices[mismatchIndexRwd]>=lvse[mismatchIndexRwd].lvsA.size() && mismatchIndexRwd<g)
			{
				mismatchIndexRwd++;
				updatedIndices[mismatchIndexRwd]++;
			}
			if (mismatchIndexRwd==g && updatedIndices[mismatchIndexRwd]>=lvse[mismatchIndexRwd].lvsA.size())
				return null;
			for (int i=0; i<mismatchIndexRwd; i++)
				updatedIndices[i]=0;
			if (mismatchIndexRwd==g)
			{
				for (int i=mismatchIndexRwd+1; i<lvse.length; i++)
					updatedIndices[i]=0;
				return updatedIndices;
			}
		}
		return updatedIndices;
	}
	
	public boolean consistentOverlaps(int [] indices)
	{
		if (indices==null)
			return false;
		//System.out.println("lvse length   "+lvse.length);
		//for (int i=0; i<lvse.length; i++)
		//	System.out.println("lvsa length "+lvse[i].lvsA.size());
		//for (int i=0; i<lvse.length; i++)
		//	System.out.println("indices "+indices[i]);
		for (int i=0; i<lvse.length-1; i++)
			if (!lvse[i].lvsA.get(indices[i]).overlaps(lvse[i+1].lvsA.get(indices[i+1])))
				return false;
		return true;
	}
	
	public GlobalVariant globalVariantReconstructor(int g)
	{
		/* old recursive code
		for (int i=0; i<lvse.length; i++)
			indices[i]=-1;
		Vector<Integer> result = new Vector();
		int o=0;
		while (o<lvse[g].lvsA.size())
		{
			result = new Vector();
			Vector<Integer> resultFwd = new Vector();
			Vector<Integer> resultRwd = new Vector();
			resultFwd.add(o);
			resultRwd.add(o);
			this.reconstruct(g,0,resultFwd,g);
			if (resultFwd.size()==(lvse.length-g))
				this.reconstructBack(g,0,resultRwd,g);
			for (int z=resultRwd.size()-1; z>0; z--)
				result.add(resultRwd.get(z));
			for (int z=0; z<resultFwd.size(); z++)
				result.add(resultFwd.get(z));
			if (result.size()==lvse.length && resultFwd.get(0)==resultRwd.get(0))
				break;
			o++;
		}
		if (result.size()!=lvse.length)
			return gv;
		for (int i=0; i<lvse.length; i++)
			indices[i]=result.get(i);
		*/
		
		//g=0;
		
		GlobalVariant gv = new GlobalVariant();
		int [] indices = new int [lvse.length];
		for (int i=0; i<indices.length; i++)
			indices[i]=0;
		System.out.print("|");
		while (!consistentOverlaps(indices))
		{
			System.out.print("\b");
			System.out.print("-");
			if (indices==null)
				return gv;
			indices=updateMatchingSet(g,indices);
			if (indices==null)
				return gv;
			System.out.print("\b");
			System.out.print("|");
		}
		
		String globalSNP = "";
		for (int i=0; i<lvse.length; i++)
		{
			if (i<(lvse.length-1))
				if (!lvse[i].lvsA.get(indices[i]).overlaps(lvse[i+1].lvsA.get(indices[i+1])))
				{
					System.out.println("\r\nglobalVariantReconstructor error (should not happen ever). Exiting program.");
					System.exit(0);
				}
			globalSNP+=lvse[i].lvsA.get(indices[i]).overlSNP1;
			globalSNP+=lvse[i].lvsA.get(indices[i]).amplSNP;
			if (i==(lvse.length-1))
				globalSNP+=lvse[i].lvsA.get(indices[i]).overlSNP2;
		}
		gv.SNP=globalSNP;
		double [] freqEst = new double [lvse.length];
		for (int i=0; i<lvse.length; i++)
		{
			freqEst[i] = lvse[i].lvsA.get(indices[i]).frequency;
		}
		LocalVariant lvGuide = lvse[g].lvsA.get(indices[g]);
		//gv.frequency=Functions.average(freqEst);
		//gv.stdevFreq=Functions.stdev(freqEst)/Math.sqrt(lvse.length);
		gv.frequency=lvGuide.frequency;
		for (int i=0; i<lvse.length; i++)
		{
			gv.frequency=Math.min(gv.frequency,lvse[i].lvsA.get(indices[i]).frequency);
		}
		for (int i=0; i<lvse.length; i++)
		{
			LocalVariant lvCurr = lvse[i].lvsA.get(indices[i]);
			lvCurr.frequency=Math.max((lvCurr.frequency-gv.frequency),0);
		}
		for (int i=0; i<lvse.length; i++)
		{
			//Vector<LocalVariant> lvA = lvse[i].lvsA;
			int j=0;
			while(j<lvse[i].lvsA.size())
			{
				if (lvse[i].lvsA.get(j).frequency<=0)
					lvse[i].lvsA.remove(j);
				else
					j++;
			}
		}
		for (int i=0; i<lvse.length; i++)
			lvse[i].sortMultinomialFrequencies();
		return gv;
	}
}
