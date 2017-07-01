import java.io.*;
import java.lang.*;
import java.util.*;
import java.math.*;

public class ClusterSet
{
	Vector<Integer> clusterSet;
	Vector<LocalVariant> localVariants;
	double length;
	String genome;
	double hErr;
	double nonhErr;
	double posteriorLikelihood;
	
	public ClusterSet()
	{
		clusterSet=null;
		localVariants=null;
		length=-1;
		genome=null;
		hErr=-1;
		nonhErr=-1;
		posteriorLikelihood=Double.MAX_VALUE;
	}
	
	public ClusterSet(Vector<Integer> c, Vector<LocalVariant> a, double l, String g, double he, double nhe)
	{
		clusterSet=c;
		localVariants=a;
		length=l;
		genome=g;
		hErr=he;
		nonhErr=nhe;
		posteriorLikelihood=Double.MAX_VALUE;
	}
	
	public void setPosterior()
	{
		posteriorLikelihood=Functions.posteriorClusteringLikelihood(clusterSet,localVariants, length, genome, hErr, nonhErr);
	}
}
		