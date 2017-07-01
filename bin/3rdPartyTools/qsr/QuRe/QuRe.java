import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class QuRe
{
	public static void main (String [] args)
	throws Exception
	{
		System.out.println("----------------------------------------------------------------------------");
		System.out.println("----------------------------------------------------------------------------");
		
		Date d = new Date();
		long starttime=d.getTime();
		
		String filename = args[0];
		if (filename.lastIndexOf('.')!=-1)
			filename=args[0].substring(0,args[0].lastIndexOf('.'));
		
		double homopolErr = 0.01d;
		double nonHomopolErr = 0.005d;
		int iterations = 3000;
		
		if (args.length>2)
		{
			homopolErr = Double.parseDouble(args[2]);
			nonHomopolErr = Double.parseDouble(args[3]);
			iterations = Integer.parseInt(args[4]);
		}
		
		int num_proc=Runtime.getRuntime().availableProcessors();
		num_proc=num_proc-1;
		num_proc=Math.max(1,num_proc);
		//num_proc=4;
		
		System.out.println("parallel processing enabled: no. of cores available = "+num_proc);
		
		ReadSet rs = new ReadSet();
		int kmer=9;
		rs.readFasta(args[0]);
		rs.readReferenceGenome(args[1],kmer);
		
		float gop=23f;
		float gep=0.3f;
		rs.setRandomScoresParallel(2000,gop,gep,num_proc);
		rs.alignParallel(num_proc,gop,gep,kmer);
		rs.setAllPvalues("BH");
		rs.removeBadReads(0.01d);
		rs.estimateBaseSet();
		rs.correctReadsParallel(0.01d, nonHomopolErr, homopolErr, num_proc);
		rs.updatePopulationStats();
		System.out.println("\t"+rs.population.size()+" reads spanning the high-coverage window");
		
		FileWriter fw = new FileWriter(filename+"_alignedReads.txt");
		//fw.write("name\tSNP_string\tmappingPosition\tstart\tstop\torientation\tsimilarity\tinsertions\tpvalue\tadjustedPvalue\r\n");
		fw.write("name\tvariations\tstart\tstop\ttrimmed_sequence\tadj_pvalue\r\n");
		for (int i=0; i<rs.population.size(); i++)
		{
			fw.write(rs.population.get(i).name+"\t");
			fw.write(rs.population.get(i).SNP_string+"\t");
			//fw.write(rs.population.get(i).mappingPosition+"\t");
			fw.write(rs.population.get(i).start+"\t");
			fw.write(rs.population.get(i).stop+"\t");
			//fw.write(rs.population.get(i).orientation+"\t");
			//fw.write(rs.population.get(i).similarity+"\t");
			//fw.write(rs.population.get(i).insertions+"\t");
			//fw.write(rs.population.get(i).pvalue+"\t");
			fw.write(rs.population.get(i).sequence+"\t");
			fw.write(rs.population.get(i).adjustedPvalue+"");
			fw.write("\r\n");
		}
		fw.close();
		
		LinkedList<Base> base_list = new LinkedList();
		Object [] key = rs.baseSet.keySet().toArray();
		Arrays.sort(key);
		for (int i=0; i<key.length; i++)
		{
			Base b = rs.baseSet.get(key[i]);
			base_list.add(b);
		}
		
		fw = new FileWriter(filename+"_snpTable.txt");
		//fw.write("position\treference\tconsensus\tA\tC\tG\tT\tdel\tcoverage\tprobcoverage\tentropy\r\n");
		fw.write("position\treference\tconsensus\tA\tC\tG\tT\tdel\tcoverage\tentropy\r\n");
		for (int i=0; i<base_list.size(); i++)
		{
			Base b = base_list.get(i);
			fw.write(b.position+"\t");
			fw.write(b.reference+"\t");
			fw.write(b.consensus+"\t");
			fw.write(b.A+"\t");
			fw.write(b.C+"\t");
			fw.write(b.G+"\t");
			fw.write(b.T+"\t");
			fw.write(b.del+"\t");
			fw.write(b.coverage+"\t");
			//boolean indel = false;
			//if (b.reference=='-' || b.consensus=='-' || Math.floor(b.position)!=b.position)
			//	indel=true;
			//fw.write(Functions.isHomopolymeric((int)Math.floor(b.position-1),rs.consensusGenomeNoIndels,indel)+"\t");
			//fw.write(b.probCoverage+"\t");
			fw.write(b.entropy+"\r\n");
			
		}
		fw.close();
		
		d = new Date();
		long stoptime1=d.getTime();
		long timePassed1=stoptime1-starttime;
		System.out.println("alignment and mapping time = "+timePassed1+" ms");
		
		System.out.println("starting Quasispecies Reconstruction (QuRe)");
		rs.estimateAmpliconsParallel(iterations, num_proc);
		
		/*
		Vector<GlobalVariant> gvFinal = new Vector();
		System.out.print("\texecuting core reconstr. algor. on "+rs.goodAmpliconSets.size()+" overl. sets ");
		for (int i=0; i<rs.goodAmpliconSets.size(); i++)
		{
			String perc = (Math.round(100f*((float)(i))/((float)(rs.goodAmpliconSets.size()))))+"% ";
			System.out.print(perc);
			rs.ampliconSet=rs.goodAmpliconSets.get(i);
			LocalVariantSetEnsemble lvse = new LocalVariantSetEnsemble(rs.population, rs.ampliconSet, rs.referenceGenome, homopolErr, nonHomopolErr, num_proc);
			//if (i==0) lvse.printToFile(filename+"_overlappingIntervalsSet.txt");
			Vector<GlobalVariant> gvFinal1 = lvse.quasispeciesReconstructor(rs.referenceGenome);
			while(gvFinal1.size()==0 && lvse.starts.length>1)
			{
				double [] newStarts = new double [rs.ampliconSet.starts.length-1];
				double [] newStops = new double [rs.ampliconSet.stops.length-1];
				if (Math.random()<0.5d)
				{
					for (int i1=0; i1<rs.ampliconSet.starts.length-1; i1++)
					{
						newStarts[i1]=rs.ampliconSet.starts[i1];
						newStops[i1]=rs.ampliconSet.stops[i1];
					}
					rs.ampliconSet.starts=newStarts;
					rs.ampliconSet.stops=newStops;
				}
				else
				{
					for (int i1=0; i1<rs.ampliconSet.starts.length-1; i1++)
					{
						newStarts[i1]=rs.ampliconSet.starts[i1+1];
						newStops[i1]=rs.ampliconSet.stops[i1+1];
					}
					rs.ampliconSet.starts=newStarts;
					rs.ampliconSet.stops=newStops;
				}
				lvse = new LocalVariantSetEnsemble(rs.population, rs.ampliconSet, rs.referenceGenome, homopolErr, nonHomopolErr, num_proc);
				lvse.printToFile(filename+"_overlappingIntervalsSet.txt");
				gvFinal1 = lvse.quasispeciesReconstructor(rs.referenceGenome);
			}
			gvFinal1 = Functions.mergeGlobalVariantSet(gvFinal1);
			Functions.setFrequencies(gvFinal1);
			for (int j=0; j<gvFinal1.size(); j++)
				gvFinal.add(gvFinal1.get(j));
			for (int o=0; o<perc.length()+8; o++)
				System.out.print("\b");
		}
		gvFinal = Functions.mergeGlobalVariantSet(gvFinal);
		Functions.setFrequencies(gvFinal);
		System.out.println("100%         ");
		System.out.println("\t\tinitial number of variants = "+gvFinal.size());
		*/
		
		System.out.print("\texecuting core reconstruction algorithm ");
		rs.ampliconSet=rs.goodAmpliconSets.get(0);
		LocalVariantSetEnsemble lvse = new LocalVariantSetEnsemble(rs.population, rs.ampliconSet, rs.referenceGenome, homopolErr, nonHomopolErr, num_proc);
		lvse.printToFile(filename+"_overlappingIntervalsSet.txt");
		Vector<GlobalVariant> gvFinal = lvse.quasispeciesReconstructor(rs.referenceGenome);
		while(gvFinal.size()==0 && lvse.starts.length>1)
		{
			double [] newStarts = new double [rs.ampliconSet.starts.length-1];
			double [] newStops = new double [rs.ampliconSet.stops.length-1];
			if (Math.random()<0.5d)
			{
				for (int i1=0; i1<rs.ampliconSet.starts.length-1; i1++)
				{
					newStarts[i1]=rs.ampliconSet.starts[i1];
					newStops[i1]=rs.ampliconSet.stops[i1];
				}
				rs.ampliconSet.starts=newStarts;
				rs.ampliconSet.stops=newStops;
			}
			else
			{
				for (int i1=0; i1<rs.ampliconSet.starts.length-1; i1++)
				{
					newStarts[i1]=rs.ampliconSet.starts[i1+1];
					newStops[i1]=rs.ampliconSet.stops[i1+1];
				}
				rs.ampliconSet.starts=newStarts;
				rs.ampliconSet.stops=newStops;
			}
			lvse = new LocalVariantSetEnsemble(rs.population, rs.ampliconSet, rs.referenceGenome, homopolErr, nonHomopolErr, num_proc);
			lvse.printToFile(filename+"_overlappingIntervalsSet.txt");
			gvFinal = lvse.quasispeciesReconstructor(rs.referenceGenome);
		}
		gvFinal= Functions.mergeGlobalVariantSet(gvFinal);
		Functions.setFrequencies(gvFinal);
		System.out.println();
		System.out.println("\t\tinitial number of variants = "+gvFinal.size());
		
		/**/
		fw = new FileWriter(filename+"_reconstructedVariantsNoClustering.txt");
		for (int i=0; i<gvFinal.size(); i++)
		{
			fw.write(">QuReNC"+i+"_"+gvFinal.get(i).frequency+"\r\n");
			fw.write(gvFinal.get(i).sequence+"\r\n");
		}
		fw.close();
		/**/
		
		System.out.print("\t\tfinal clustering (random search + BIC selection) ");
		if (gvFinal.size()>1)
			gvFinal = Functions.correct(gvFinal, rs.ampliconSet.starts[0], rs.ampliconSet.stops[rs.ampliconSet.stops.length-1], rs.referenceGenome, homopolErr, nonHomopolErr, iterations, num_proc);
		gvFinal = Functions.mergeGlobalVariantSet(gvFinal);
		Functions.setFrequencies(gvFinal);
		System.out.println("\t\tfinal number of variants = "+gvFinal.size());
		
		fw = new FileWriter(filename+"_reconstructedVariants.txt");
		for (int i=0; i<gvFinal.size(); i++)
		{
			fw.write(">QuRe"+i+"_"+gvFinal.get(i).frequency+"\r\n");
			fw.write(gvFinal.get(i).sequence+"\r\n");
		}
		fw.close();
		
		d = new Date();
		long stoptime2=d.getTime();
		long timePassed2=stoptime2-stoptime1;
		long timePassed3=stoptime2-starttime;
		System.out.println("amplicon estimation and quasispecies reconstruction time = "+timePassed2+" ms");
		System.out.println("total time employed = "+timePassed3+" ms");
		
		System.out.println("----------------------------------------------------------------------------");
		System.out.println("----------------------------------------------------------------------------");
	}
}
