import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class AmpliconSet
{
	double [] starts;
	double [] stops;
	double minReadCoverage;
	double overallReadCoverage;
	double minOverlapDiversity;
	double overallOverlapDiversity;
	double numNonZeroOverlapDiversity;
	double numOverlaps;
	double minOverlapLength;
	double overallOverlapLength;
	double minAmpliconLength;
	double overallAmpliconLength;
	double minReadCoverageProbability;
	double overallReadCoverageProbability;
	double minOverlapDiversityProbability;
	double overallOverlapDiversityProbability;
	double numNonZeroOverlapDiversityProbability;
	double numOverlapsProbability;
	double minOverlapLengthProbability;
	double overallOverlapLengthProbability;
	double minAmpliconLengthProbability;
	double overallAmpliconLengthProbability;
	double score;
	
	public AmpliconSet()
	{
		starts=null;
		stops=null;
		minReadCoverage=0;
		overallReadCoverage=0;
		minOverlapDiversity=0;
		overallOverlapDiversity=0;
		numNonZeroOverlapDiversity=0;
		numOverlaps=0;
		minOverlapLength=0;
		overallOverlapLength=0;
		minAmpliconLength=0;
		overallAmpliconLength=0;
		minReadCoverageProbability=0;
		overallReadCoverageProbability=0;
		minOverlapDiversityProbability=0;
		overallOverlapDiversityProbability=0;
		numNonZeroOverlapDiversityProbability=0;
		numOverlapsProbability=0;
		minOverlapLengthProbability=0;
		overallOverlapLengthProbability=0;
		minAmpliconLengthProbability=0;
		overallAmpliconLengthProbability=0;
		score=-1E100;
	}
	
	public AmpliconSet(double [] sta, double [] sto, Vector<Read> population, double alignStart, double alignStop)
	{
		score=0;
		starts=sta;
		stops=sto;
		Vector<Read> [] covering = new Vector [starts.length];
		double mrc = 0;
		double [] rcv = new double [starts.length];
		for (int j=0; j<starts.length; j++)
		{
			covering[j]=new Vector();
			double rc=0;
			for (int i=0; i<population.size(); i++)
			{
				Read r = population.get(i);
				if (r.spans(starts[j],stops[j]))
				{
					rc++;
					covering[j].add(r);
				}
			}
			if (j==0)
				mrc=rc;
			if (rc<mrc)
				mrc=rc;
			rcv[j]=rc;
		}
		minReadCoverage=mrc;
		overallReadCoverage=Functions.average(rcv);
		double modp = 0;
		double [] modp2 = new double[starts.length-1];
		for (int j=1; j<starts.length; j++)
		{
			/*
			Hashtable hh = new Hashtable();
			double cover=1;
			for (int i=0; i<population.size(); i++)
			{
				Read r = population.get(i);
				if (r.spans(starts[j],stops[j-1]))
				{
					cover++;
					String s = r.getSNPString(starts[j],stops[j-1]);
					if (hh.get(s)==null)
						hh.put(s,1);
				}
			}
			double odp=(double)(hh.keySet().toArray().length+1)/cover;
			*/
			Vector<Read> allOverlappingReads = new Vector();
			for (int i=0; i<covering[j-1].size(); i++)
				allOverlappingReads.add(covering[j-1].get(i));
			for (int i=0; i<covering[j].size(); i++)
				allOverlappingReads.add(covering[j].get(i));
			int subSample = 3333;
			double odp=0;
			if (minReadCoverage>0 && allOverlappingReads.size()>0)
			{
				if ((double)allOverlappingReads.size()*((double)allOverlappingReads.size()-1d)/2d<subSample)
				{
					for (int i=0; i<allOverlappingReads.size(); i++)
						for (int k=i+1; k<allOverlappingReads.size(); k++)
						{
							Read r1 = allOverlappingReads.get(i);
							Read r2 = allOverlappingReads.get(k);
							double dista = r1.distance(r2,starts[j],stops[j-1])/(stops[j-1]-starts[j]);
							odp+=dista;
						}
					odp=odp/((double)allOverlappingReads.size()*((double)allOverlappingReads.size()-1d)/2d);
					//System.out.println("full "+odp);
				}
				else
				{
					for (int i=0; i<subSample; i++)
					{
						Read r1 = allOverlappingReads.get((int)(Math.random()*(allOverlappingReads.size()-1)));
						Read r2 = allOverlappingReads.get((int)(Math.random()*(allOverlappingReads.size()-1)));
						double dista = r1.distance(r2,starts[j],stops[j-1])/(stops[j-1]-starts[j]);
						odp+=dista;
					}
					odp=odp/(double)subSample;
					//System.out.println("subSample "+odp);
				}
			}
			if (j==1)
			{
				modp=odp;
			}
			else
			{
				if (odp<modp)
					modp=odp;
			}
			modp2[j-1]=odp;
		}
		minOverlapDiversity=modp;
		overallOverlapDiversity=Functions.average(modp2);
		double [] ol = new double [modp2.length];
		for (int j=0; j<modp2.length; j++)
		{
			ol[j]=stops[j]-starts[j+1];
			if (j==0)
				minOverlapLength=stops[j]-starts[j+1];
			if ((stops[j]-starts[j+1])<minOverlapLength)
				minOverlapLength=stops[j]-starts[j+1];
			if (modp2[j]>0)
				numNonZeroOverlapDiversity++;
		}
		numNonZeroOverlapDiversity=numNonZeroOverlapDiversity/(double)starts.length;
		numOverlaps=starts.length;
		overallOverlapLength=Functions.average(ol);
		double [] al = new double [starts.length];
		for (int j=0; j<starts.length; j++)
		{
			al[j]=stops[j]-starts[j];
		}
		minAmpliconLength=Functions.min(al);
		overallAmpliconLength=Functions.average(al);
	}
	
	public void setStats(double [] avgReadCover, double [] avgOverlapDiversity, double [] avgOverallOverlapDiversity, double [] avgNumNonZeroOverlapDiversity, double [] avgOverallReadCoverage, double [] avgNumOverlaps, double [] avgMinOverlapLength, double [] avgOverallOverlapLength, double [] avgMinAmpliconLength, double [] avgOverallAmpliconLength)
	{
		minReadCoverageProbability=0;
		overallReadCoverageProbability=0;
		minOverlapDiversityProbability=0;
		overallOverlapDiversityProbability=0;
		numNonZeroOverlapDiversityProbability=0;
		numOverlapsProbability=0;
		minOverlapLengthProbability=0;
		overallOverlapLengthProbability=0;
		for (int i=0; i<avgReadCover.length; i++)
		{
			if (avgReadCover[i]<=minReadCoverage) minReadCoverageProbability++;
			if (avgOverallReadCoverage[i]<=overallReadCoverage) overallReadCoverageProbability++;
			if (avgOverlapDiversity[i]<=minOverlapDiversity) minOverlapDiversityProbability++;
			if (avgOverallOverlapDiversity[i]<=overallOverlapDiversity) overallOverlapDiversityProbability++;
			if (avgNumNonZeroOverlapDiversity[i]<=numNonZeroOverlapDiversity) numNonZeroOverlapDiversityProbability++;
			if (avgNumOverlaps[i]>=numOverlaps) numOverlapsProbability++;
			if (avgMinOverlapLength[i]<=minOverlapLength) minOverlapLengthProbability++;
			if (avgOverallOverlapLength[i]<=overallOverlapLength) overallOverlapLengthProbability++;
			if (avgMinAmpliconLength[i]<=minAmpliconLength) minAmpliconLengthProbability++;
			if (avgOverallAmpliconLength[i]<=overallAmpliconLength) overallAmpliconLengthProbability++;
		}
		minReadCoverageProbability=minReadCoverageProbability/avgReadCover.length;
		overallReadCoverageProbability=overallReadCoverageProbability/avgReadCover.length;
		minOverlapDiversityProbability=minOverlapDiversityProbability/avgReadCover.length;
		overallOverlapDiversityProbability=overallOverlapDiversityProbability/avgReadCover.length;
		numNonZeroOverlapDiversityProbability=numNonZeroOverlapDiversityProbability/avgReadCover.length;
		numOverlapsProbability=numOverlapsProbability/avgReadCover.length;
		minOverlapLengthProbability=minOverlapLengthProbability/avgReadCover.length;
		overallOverlapLengthProbability=overallOverlapLengthProbability/avgReadCover.length;
		minAmpliconLengthProbability=minAmpliconLengthProbability/avgReadCover.length;
		overallAmpliconLengthProbability=overallAmpliconLengthProbability/avgReadCover.length;
		/*
		minReadCoverageProbability=Functions.zetaStandardProbability((minReadCoverage - Functions.average(avgReadCover))/Functions.stdev(avgReadCover));
		overallReadCoverageProbability=Functions.zetaStandardProbability((overallReadCoverage - Functions.average(avgOverallReadCoverage))/Functions.stdev(avgOverallReadCoverage));
		minOverlapDiversityProbability=Functions.zetaStandardProbability((minOverlapDiversity - Functions.average(avgOverlapDiversity))/Functions.stdev(avgOverlapDiversity));
		overallOverlapDiversityProbability=Functions.zetaStandardProbability((overallOverlapDiversity - Functions.average(avgOverallOverlapDiversity))/Functions.stdev(avgOverallOverlapDiversity));
		numNonZeroOverlapDiversityProbability=Functions.zetaStandardProbability((numNonZeroOverlapDiversity - Functions.average(avgNumNonZeroOverlapDiversity))/Functions.stdev(avgNumNonZeroOverlapDiversity));
		*/
		score=0;
		score=score+Math.log(Math.max(1E-100,minReadCoverageProbability));
		score=score+Math.log(Math.max(1E-100,overallReadCoverageProbability));
		score=score+Math.log(Math.max(1E-100,minOverlapDiversityProbability));
		score=score+Math.log(Math.max(1E-100,overallOverlapDiversityProbability));
		score=score+Math.log(Math.max(1E-100,numNonZeroOverlapDiversityProbability));
		score=score+Math.log(Math.max(1E-100,numOverlapsProbability));
		score=score+Math.log(Math.max(1E-100,minOverlapLengthProbability));
		score=score+Math.log(Math.max(1E-100,overallOverlapLengthProbability));
		score=score+Math.log(Math.max(1E-100,minAmpliconLengthProbability));
		score=score+Math.log(Math.max(1E-100,overallAmpliconLengthProbability));
		if (Double.isNaN(score) || Double.isInfinite(score) || Double.isInfinite(-score))
			score=-1E100;
	}
	
	public boolean isBetter(AmpliconSet b)
	{
		return (this.score>b.score);
	}
	
	public boolean checkConsistency()
	{
		if (starts.length==1)
			return true;
		for (int i=0; i<(starts.length-1); i++)
		{
			double start1=starts[i];
			double start2=starts[i+1];
			double stop1=stops[i];
			double stop2=stops[i+1];
			if (!(start1<start2 && stop1<stop2 && start1<stop1 && start2<stop2 && start2<stop1))
				return false;
		}
		return true;
	}
	
	public boolean checkNonMutualOverlapsConsistency()
	{
		if (starts.length==1)
			return true;
		if (starts.length==2)
		{
			double start1=starts[0];
			double start2=starts[1];
			double stop1=stops[0];
			double stop2=stops[1];
			return ( (start2<stop1) && (start1<start2) && (stop1<stop2) );
		}
		for (int i=0; i<(starts.length-2); i++)
		{
			double start1=starts[i];
			double start2=starts[i+1];
			double start3=starts[i+2];
			double stop1=stops[i];
			double stop2=stops[i+1];
			double stop3=stops[i+2];
			if 
			(!(
				(stop1<start3) &&
				(start1<start2) &&
				(start2<start3) &&
				(stop1<stop2) &&
				(stop2<stop3) &&
				(start2<stop1) &&
				(start3<stop2)
			))
			return false;
		}
		return true;
	}
	
	public void initialiseRandomAmpliconSet(double avgReadLength, double stdReadLength, double alignStart, double alignStop, Vector<Read> pop)
	{
		double[] starts1;
		double[] stops1;
		double ch=Math.random();
		if (ch<0.9)
		{
			LinkedList<Double> startL = new LinkedList();
			LinkedList<Double> stopL = new LinkedList();
			startL.add(alignStart);
			Random rndA = new Random();
			double rngA = rndA.nextGaussian()*stdReadLength+(avgReadLength-3*stdReadLength);
			if (stdReadLength<5)
				{
					rngA = rndA.nextGaussian()*5+(avgReadLength-3*5);
					if (Math.random()<0.1) {rngA = rndA.nextGaussian()*Math.sqrt(avgReadLength)+avgReadLength/2;}
				}
			rngA = rndA.nextGaussian()*Math.sqrt(avgReadLength)+avgReadLength/1.5d;
			if (Math.random()<0.1) {rngA = Math.random()*avgReadLength;}
			if (rngA<50)
				rngA = Math.max(50,Math.random()*avgReadLength);
			stopL.add(alignStart+rngA);
			while (true)
			{
				if (stopL.getLast().doubleValue()==alignStop)
					break;
				rngA = rndA.nextGaussian()*stdReadLength+(avgReadLength-3*stdReadLength);
				if (stdReadLength<5)
				{
					rngA = rndA.nextGaussian()*5+(avgReadLength-3*5);
					if (Math.random()<0.1) {rngA = rndA.nextGaussian()*Math.sqrt(avgReadLength)+avgReadLength/2;}
				}
				if (Math.random()<0.1) {rngA = Math.random()*avgReadLength;}
				if (rngA<50)
					rngA = Math.max(50,Math.random()*avgReadLength);
				double minStart = Math.max(startL.getLast(),stopL.getLast()-rngA+1);
				double actStart = minStart+Math.random()*(stopL.getLast()-minStart);
				double actStop = Math.max(stopL.getLast()+1,actStart+rngA);
				actStop = Math.min(alignStop,actStop);
				startL.add(actStart);
				stopL.add(actStop);
			}
			starts1 = new double[startL.size()];
			stops1 = new double[startL.size()];
			for (int l=0; l<starts1.length; l++)
			{
				starts1[l]=startL.get(l);
				stops1[l]=stopL.get(l);
			}
		}
		else
		{
			Random rndA = new Random();
			double rngA = rndA.nextGaussian()*stdReadLength+(avgReadLength-3*stdReadLength);
			if (stdReadLength<5)
			{
				rngA = rndA.nextGaussian()*5+(avgReadLength-3*5);
				if (Math.random()<0.1) {rngA = rndA.nextGaussian()*Math.sqrt(avgReadLength)+avgReadLength/2;}
			}
			if (Math.random()<0.1) {rngA = Math.random()*avgReadLength;}
			if (rngA<50)
				rngA = Math.max(50,Math.random()*avgReadLength);
			double window = rngA;
			if (alignStart+window>alignStop)
				window=(alignStop-alignStart)/3;
			double step = Math.max(5,window*Math.random());
			LinkedList<Double> startL = new LinkedList();
			LinkedList<Double> stopL = new LinkedList();
			startL.add(alignStart);
			stopL.add(alignStart+window+step);
			while(stopL.getLast().doubleValue()<alignStop)
			{
				startL.add(startL.getLast()+step);
				stopL.add(stopL.getLast()+step);
				if (stopL.getLast().doubleValue()>alignStop)
					{
						stopL.removeLast();
						stopL.add(alignStop);
					}
			}
			starts1 = new double[startL.size()];
			stops1 = new double[startL.size()];
			for (int l=0; l<starts1.length; l++)
			{
				starts1[l]=startL.get(l);
				stops1[l]=stopL.get(l);
			}
		}
		AmpliconSet amplicon = new AmpliconSet(starts1,stops1,pop,alignStart,alignStop);
		starts=amplicon.starts;
		stops=amplicon.stops;
		minReadCoverage=amplicon.minReadCoverage;
		overallReadCoverage=amplicon.overallReadCoverage;
		minOverlapDiversity=amplicon.minOverlapDiversity;
		overallOverlapDiversity=amplicon.overallOverlapDiversity;
		numNonZeroOverlapDiversity=amplicon.numNonZeroOverlapDiversity;
		numOverlaps=amplicon.numOverlaps;
		minOverlapLength=amplicon.minOverlapLength;
		overallOverlapLength=amplicon.overallOverlapLength;
		minAmpliconLength=amplicon.minAmpliconLength;
		overallAmpliconLength=amplicon.overallAmpliconLength;
		minReadCoverageProbability=amplicon.minReadCoverageProbability;
		overallReadCoverageProbability=amplicon.overallReadCoverageProbability;
		minOverlapDiversityProbability=amplicon.minOverlapDiversityProbability;
		overallOverlapDiversityProbability=amplicon.overallOverlapDiversityProbability;
		numNonZeroOverlapDiversityProbability=amplicon.numNonZeroOverlapDiversityProbability;
		minOverlapLengthProbability=amplicon.minOverlapLengthProbability;
		overallOverlapLengthProbability=amplicon.overallOverlapLengthProbability;
		score=amplicon.score;
	}
}
	