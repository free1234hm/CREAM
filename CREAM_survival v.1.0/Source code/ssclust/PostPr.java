package ssclust;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

public class PostPr {
	public double post(double[] express, double[] mean, double zeta, double varht){
		double temp=0;
		int m = express.length;

		double tempvar = Math.pow(10, zeta)*varht;
		double rho=tempvar/(tempvar+varht);
		double[][] sigma = new double[m][m];
		for(int i=0;i<m;i++){
			for(int j=0;j<m;j++) sigma[i][j] = rho;
		}
		for(int i=0;i<m;i++){
			sigma[i][i] = sigma[i][i]+(1-rho);
		}
		for(int i=0;i<m;i++){
			for(int j=0;j<m;j++) sigma[i][j] = sigma[i][j]*(tempvar+varht);
		}
		temp = dmvnorm(express, mean, sigma); 
		/**xx is gene expression level,
		 * mean is mean curve
		 * sigma is covariance 
		 */
		return temp;
	}
	
	public double dmvnorm(double[] xx1, double[] mean, double[][] sigma1){
		double[] xx = new double[xx1.length];
		for(int i=0;i<xx.length;i++) xx[i] = xx1[i];
		double[][] sigma = new double[sigma1.length][sigma1[0].length];
		for(int i=0;i<sigma.length;i++){
			for(int j=0;j<sigma[0].length;j++){
				sigma[i][j] = sigma1[i][j];
			}
		}
		int p = xx.length;
		double logretval=0;
		/*²»¿¼ÂÇxxÈ±Ê§Öµ
		List<Integer> unmiss = new ArrayList<Integer>(); 
		for(int i=0;i<xx.length;i++){
			if(xx[i] != null) unmiss.add(i);
		}
		int p = unmiss.size();
		if(unmiss.size() < xx.length){
			Double[] xx1 = new Double[p];
			double[] mean1 = new double[p];
			double[][] sigma1 = new double[p][p];
			for(int i=0;i<p;i++){
				xx1[i] = xx[unmiss.get(i)];
				mean1[i] = mean[unmiss.get(i)];
				for(int j=0;j<p;j++){
					sigma1[i][j] = sigma[unmiss.get(i)][unmiss.get(j)];
				}
			}
			xx = xx1;
			mean = mean1;
			sigma = sigma1;
		}
		*/	
		SpineCurve sp = new SpineCurve();
		xx = sp.getspline(xx, mean);
		boolean error = false;
		if(p != sigma[0].length){
			error = true;
		}
		double[] work = new double[p];
		Integer[] jpvt = new Integer[p];
		int info = Cholesky.a(sigma, sigma.length, sigma[0].length, work, jpvt, 0);
		for(int i=0;i<sigma.length;i++){
			if(sigma[i][i]<=0) sigma[i][i] = 2.220446e-16;
		}
		
		if(error){
			int sum = 0;
			for(int i=0;i<p;i++){
				if(xx[i] == mean[i]) sum++;
			}
			if(sum == p) {
				logretval = 0;
			}else{
				logretval =  Double.NEGATIVE_INFINITY;
			}
			
		}else{
			double[] deviation = new double[p];
			for(int i=0;i<p;i++) deviation[i] = xx[i] - mean[i];
			Solve solve = new Solve();
			double[] tmp = solve.backsolve(sigma, deviation, true);
	        double rss = 0;
	        for(int i=0;i<p;i++) rss = rss+tmp[i]*tmp[i];
	        for(int i=0;i<p;i++) logretval = logretval - Math.log(sigma[i][i]);
	        logretval = logretval - 0.5 * p * Math.log(2 * Math.PI) - 0.5*rss;    
		}
		//logretval = Math.pow(Math.E, logretval);
		return logretval;
	}

	/*
	public static void main(String[] args) {
		double[] x = {2.0, 3.3, 3.0, 5.4, 4.6};
		double[] mean = {12.0, 13.2, 13.3, 15.0, 14.5};
		double[][] cov = {{1.5, 1.2, 1.2, 1.2, 1.2},
				          {1.2, 1.5, 1.2, 1.2, 1.2},
				          {1.2, 1.2, 1.5, 1.2, 1.2},
				          {1.2, 1.2, 1.2, 1.5, 1.2},
				          {1.2, 1.2, 1.2, 1.2, 1.5}};
		//double pr = post(x,mean,zeta, varht);
		
		double pr = dmvnorm(x,mean,cov);
		
	}
	*/

}
