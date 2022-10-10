package bsm.core;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

import JSci.maths.statistics.FDistribution;

public class Anova {

	Variance variance = new Variance();
	Mean meanUtil = new Mean();
	Sum sumUtil = new Sum();

	private List<double[]> list=new ArrayList<double[]>();

	public void addE(double[] e) {
		list.add(e);
	}

	public List<Double> getVariance() {
		List<Double> varianceMap = new ArrayList<Double>();
		for (int i=0;i<list.size();i++){
			double[] e = list.get(i);
			varianceMap.add(variance.evaluate(e));
		}
		return varianceMap;
	}

	public List<Double> getMean() {
		List<Double> meanMap = new ArrayList<Double>();
		for (int i=0;i<list.size();i++){
			double[] e = list.get(i);
			meanMap.add(meanUtil.evaluate(e));
		}
		return meanMap;
	}

	public double getSumMean() {
		double sum = 0;
		int n = 0;
		for (int i=0;i<list.size();i++){
			double[] e = list.get(i);
			sum = sum + sumUtil.evaluate(e);
			n = n + e.length;
		}
		return sum / n;
	}
	
	public int getSumNum() {
		int n = 0;
		for (int i=0;i<list.size();i++){
			double[] e = list.get(i);
			n = n + e.length;
		}
		return  n;
	}


	public double getMSTR() {
		double sumMean = getSumMean();
		double numerator=0;
		for (int i=0;i<list.size();i++){
			double[] e = list.get(i);
			double mean=meanUtil.evaluate(e);
			numerator = numerator+e.length*(mean-sumMean)*(mean-sumMean);
		}
		return numerator/(list.size()-1);
	}

	public double getMSE() {
		double numerator=0;
		for (int i=0;i<list.size();i++){
			double[] e = list.get(i);
			double v = variance.evaluate(e);
			numerator = numerator+(e.length-1)*v;
		}
		double denominator=getSumNum()-list.size();
		return numerator/denominator;
	}
	
	public double PValue() {
		double MSTR= getMSTR();
		double MSE=getMSE();
		double f=MSTR/MSE;
		double free1=list.size()-1;
		double free2=getSumNum()-list.size();
		FDistribution fd=new FDistribution(free1, free2);
		return 1-fd.cumulative(f);
	}
}

