package bsm;

import org.apache.commons.math3.stat.inference.TestUtils;
import org.apache.commons.math3.stat.inference.MannWhitneyUTest;

import bsm.core.Anova;
import smile.stat.distribution.GaussianDistribution;
import smile.stat.hypothesis.FTest;
import smile.stat.hypothesis.KSTest;

public class GetDEList {
	public double[][] test(Double[][] controldata, Double[][] casedata, int stat) throws Exception {
	    
		double[][] delist = new double[controldata.length][3];
	    
	    if(stat == 0){
	    	double[] controlvalue = new double[controldata[0].length];
	 	    double[] casevalue = new double[casedata[0].length];
	    	for(int i=0;i<delist.length;i++){
	    		for(int m=0;m<controlvalue.length;m++) controlvalue[m] = controldata[i][m];
	    		for(int m=0;m<casevalue.length;m++) casevalue[m] = casedata[i][m];
				double mean1 = 0;
				double mean2 = 0;
				for(double v:controlvalue) mean1 += v/controlvalue.length;
				for(double v:casevalue) mean2 += v/casevalue.length;
				double fc = mean2-mean1;
				
				Anova anovo=new Anova();
				anovo.addE(controlvalue);
				anovo.addE(casevalue);
				
				delist[i][0] = anovo.PValue();
				delist[i][2] = fc;
			}
	    	for(int i=0;i<delist.length;i++){
	    		int rank = getrank(delist[i][0], delist);
	    		delist[i][1] = (delist[i][0]*delist.length)/rank;
	    	}
	    	
	    }else if(stat == 1){
	    	
	    	for(int i=0;i<delist.length;i++){
	    		double[] controlvalue = new double[controldata[0].length];
		 	    double[] casevalue = new double[casedata[0].length];
	    		
	    		for(int m=0;m<controlvalue.length;m++) controlvalue[m] = controldata[i][m];
	    		for(int m=0;m<casevalue.length;m++) casevalue[m] = casedata[i][m];
				double mean1 = 0;
				double mean2 = 0;
				for(double v:controlvalue) mean1 += v/controlvalue.length;
				for(double v:casevalue) mean2 += v/casevalue.length;
				double fc = mean2-mean1;
				delist[i][0] = TestUtils.tTest(controlvalue, casevalue);
				delist[i][2] = fc;
			}
	    	for(int i=0;i<delist.length;i++){
	    		int rank = getrank(delist[i][0], delist);
	    		delist[i][1] = (delist[i][0]*delist.length)/rank;
	    	}
	    	
	    }else if(stat == 2){
	    	MannWhitneyUTest sumrank = new MannWhitneyUTest();
	    	double[] controlvalue = new double[controldata[0].length];
	 	    double[] casevalue = new double[casedata[0].length];
	    	
	    	for(int i=0;i<delist.length;i++){
	    		for(int m=0;m<controlvalue.length;m++) controlvalue[m] = controldata[i][m];
	    		for(int m=0;m<casevalue.length;m++) casevalue[m] = casedata[i][m];
				double mean1 = 0;
				double mean2 = 0;
				for(double v:controlvalue) mean1 += v/controlvalue.length;
				for(double v:casevalue) mean2 += v/casevalue.length;
				double fc = mean2-mean1;
				delist[i][0] = sumrank.mannWhitneyUTest(controlvalue, casevalue);
				delist[i][2] = fc;
			}
	    	for(int i=0;i<delist.length;i++){
	    		int rank = getrank(delist[i][0], delist);
	    		delist[i][1] = (delist[i][0]*delist.length)/rank;
	    	}
	    }else{
	    	MannWhitneyUTest sumrank = new MannWhitneyUTest();
	    	double[] controlvalue = new double[controldata[0].length];
	 	    double[] casevalue = new double[casedata[0].length];
	 	    int tcount = 0;
			int ucount = 0;
	    	for(int i=0;i<delist.length;i++){
	    		for(int m=0;m<controlvalue.length;m++) controlvalue[m] = controldata[i][m];
	    		for(int m=0;m<casevalue.length;m++) casevalue[m] = casedata[i][m];
				double mean1 = 0;
				double mean2 = 0;
				for(double v:controlvalue) mean1 += v/controlvalue.length;
				for(double v:casevalue) mean2 += v/casevalue.length;
				double fc = mean2-mean1;
				GaussianDistribution gaussian1 = new GaussianDistribution(controlvalue);
				KSTest ks1 = KSTest.test(controlvalue, gaussian1);
				double p1 = ks1.pvalue;
				
				GaussianDistribution gaussian2 = new GaussianDistribution(casevalue);
				KSTest ks2 = KSTest.test(casevalue, gaussian2);
				double p2 = ks2.pvalue;
				
				if(p1>0.05 && p2>0.05){ //Normality test
					FTest ff = FTest.test(controlvalue, casevalue);
					double fpvalue = ff.pvalue;
					if(fpvalue > 0.05){  //Homogeneity of variance test
						delist[i][0] = TestUtils.tTest(controlvalue, casevalue);
						delist[i][2] = fc;
						tcount++;
					}else{
						delist[i][0] = sumrank.mannWhitneyUTest(controlvalue, casevalue);
						delist[i][2] = fc;
						ucount++;
					}
				}else{
					delist[i][0] = sumrank.mannWhitneyUTest(controlvalue, casevalue);
					delist[i][2] = fc;
					ucount++;
				}
			}
	    	System.out.println("ttest: "+tcount+ " utest: "+ucount);
	    	for(int i=0;i<delist.length;i++){
	    		int rank = getrank(delist[i][0], delist);
	    		delist[i][1] = (delist[i][0]*delist.length)/rank;
	    	}
	    }
	    return delist;
	}
	
	private int getrank(double value, double[][] list){
		int rank = 1;
		for(int i=0;i<list.length;i++){
			if(value > list[i][0]){
				rank++;
			}
		}
		return rank;
	}
	
	/*
	 * F-test检验两组数的方差是否差异（p<0.05）；anova使用了f检验，但检验两组数的均值是否差异。
	public static void main(String[] args) throws Exception {
		 double[] x = new double[10];
		    for (int i = 0; i < x.length; i++) {
		        x[i] = Math.random()*100+100;
		    }
		    double[] y = new double[10];
		    for (int i = 0; i < y.length; i++) {
		        y[i] = Math.random()*100+30;
		    }
		    FTest ff = FTest.test(x, y);
			double fpvalue = ff.pvalue;
		    System.out.println(fpvalue);
		    
		    Anova anovo=new Anova();
			anovo.addE(x);
			anovo.addE(y);
			System.out.println(anovo.PValue());
	}
	*/
}
