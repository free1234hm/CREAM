package ssclust;

import java.util.ArrayList;
import java.util.List;

public class SpineCurve {

	public double[] getspline(double[] mydata, double[] mean) { 
		int[] yy = new int[mean.length];
		List<Double> a = new ArrayList<Double>();
		List<Double> b = new ArrayList<Double>();
		for(int i=0;i<yy.length;i++) {
			yy[i] = i+1;
			a.add((double)(i+1));
			b.add(mydata[i]);
		}
		double x[] = new double[a.size()];
		double y[] = new double[a.size()];
		for(int i=0;i<x.length;i++){
			x[i] = a.get(i);
			y[i] = b.get(i);
		}
		int n=x.length;
		double x1[]=new double[n+2];
		double y1[]=new double[n+2];
		x1[0]=y1[0]=0;
		x1[n+1]=y1[n+1]=0;
		for(int i=1;i<n+1;i++)
			{ 
			   x1[i]=x[i-1];
			   y1[i]=y[i-1];
			}
		Spline ca=new Spline(x1 ,y1, n+2);
		double s[]=ca.s(yy);
		double[] result = new double[yy.length];
		for (int i=0;i<yy.length;i++)
		{
		result[i] = s[yy[i]];	
		}
		return result;
	}
	
    /*
	public static void main(String []args)
	{
	Double y[]={1.0, 4.0, null, 10.5,24.7,64.0,null};
	double[] mean = {1,2,3,4,5,6,7};
	
	Double s[]=getspline(y, mean);
	for (int i=0;i<s.length;i++)
	{
	System.out.println(s[i]);	
	}
	}
	*/
	
	}
