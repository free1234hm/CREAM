package ssclust;

import java.util.ArrayList;
import java.util.List;

public class Gettrc {
	double mchpr = 2.220446e-16;
	public double gettrc(double[][] myweight) {
		
		Integer env = myweight.length;
		Integer time = myweight[0].length;
		double[] weight;
		Integer[] tm;
		List<Integer> tim = new ArrayList<Integer>();
		List<Double> weightt = new ArrayList<Double>();
		
			for(int i = 0;i<env;i++){ 
				double sum=0;
				for(int j=0;j<time;j++){ //除了初始节点其余weight均为0时跳过此基因
					sum += myweight[i][j];
				}
				if(sum == 0) continue;
				for(int j=0;j<time;j++){					
					if(myweight[i][j] > 0){
						weightt.add(myweight[i][j]);
						tim.add(j+1);
					}
				}
			}

			weight = new double[weightt.size()];
			tm = new Integer[tim.size()];
			for(int i=0;i<weight.length;i++){
				weight[i] = weightt.get(i);
				tm[i] = tim.get(i);
			}
			Integer[] id_basis = new Integer[time];	
			for(int i=0;i<time;i++) id_basis[i] = i+1;
		
		
		double env_max=tm[0],env_min=tm[0];
		for(int i=0;i<tm.length;i++){
			if(env_max<tm[i]) env_max = tm[i];
			if(env_min>tm[i]) env_min = tm[i];
		}
		double distance = env_max - env_min;
		env_max = env_max + distance*0.3;
		env_min = env_min - distance*0.3;
		double[] x_basis = new double[id_basis.length];
		for(int i=0;i<x_basis.length;i++){
			x_basis[i] = tm[id_basis[i]-1];
		}
		double[][] s = mkphi(tm, env_min, env_max);
		double[][] r = mkrk(tm, x_basis, env_min, env_max);
		double[][] q = new double[id_basis.length][r[0].length];
		for(int i=0;i<q.length;i++){
			for(int j=0;j<q[0].length;j++){
				q[i][j] = r[id_basis[i]-1][j];
			}	
		}
		if(weight != null && weight.length>0){ 
			double[] weight1 = new double[weight.length];
			for(int i=0;i<weight1.length;i++) weight1[i] = weight[i];
			for(int i=0;i<s.length;i++){
				for(int j=0;j<s[0].length;j++) s[i][j] = s[i][j] * weight1[i];
			}
			for(int i=0;i<r.length;i++){
				for(int j=0;j<r[0].length;j++) r[i][j] = r[i][j] * weight1[i];
			}
		}
		
		int nobs = r.length;  //基因数*时间数
		int nxi = r[0].length; //时间数
		int nnull = s[0].length;
		int nxiz = nxi;   //时间数
		int nn = nxiz + nnull;  //时间数+2
		
		double tmp = 0;
		double theta = 0;
		for(int i=0;i<r.length;i++){
			for(int j=0;j<r[0].length;j++){
				tmp = tmp + r[i][j]*r[i][j];
			}
		}
		if(s == null || s.length == 0){
			theta = 0;
		}else{
			double sum_s = 0;
			for(int i=0;i<s.length;i++){
				for(int j=0;j<s[0].length;j++){
				sum_s = sum_s + s[i][j]*s[i][j];
				}
				}
			theta = Math.log10(sum_s/nnull/tmp*nxi)/2;
		}
		
		double[][] q_wk;
			q_wk = new double[nxi][nxi];
			for(int i=0;i<nxi;i++){
				for(int j=0;j<nxi;j++){
					q_wk[i][j] = Math.pow(10, theta)*q[i][j];
				}
			}
			double[][] tq_wk = transpose(q_wk);
			for(int i=0;i<nxiz;i++){
				for(int j=0;j<nxiz;j++){
					q_wk[i][j] = (q_wk[i][j] + tq_wk[i][j])/2;
				}
			}
		
		
		double[][] qq_wk = new double[q_wk.length][q_wk[0].length];
		for(int i=0;i<q_wk.length;i++){
			for(int j=0;j<q_wk[0].length;j++) qq_wk[i][j] = q_wk[i][j];
		}
		for(int i=0;i<q_wk.length;i++) qq_wk[i][i] = qq_wk[i][i]+mchpr;
		double[] work = new double[qq_wk.length];
		Integer[] jpvt = new Integer[qq_wk.length];
		for(int i=1;i<=jpvt.length;i++) jpvt[i-1] = 0;
		int rank = Cholesky.a(qq_wk, qq_wk.length, qq_wk[0].length, work, jpvt, 1);
		double[][] pivot = new double[r.length][r[0].length];
		for(int i=0;i<r[0].length;i++){
			for(int j=0;j<r.length;j++) pivot[j][i] = Math.pow(10, theta)*r[j][jpvt[i]-1];
		}
		double[][] sr1 = new double[s.length][s[0].length+pivot[0].length];
		for(int i=0;i<s[0].length;i++){
			for(int j=0;j<s.length;j++) sr1[j][i] = s[j][i];
		}
		for(int i=0;i<pivot[0].length;i++){
			for(int j=0;j<s.length;j++) sr1[j][i+s[0].length] = pivot[j][i];
		}
		
		double[][] q_wk2 = new double[nxiz][nnull+qq_wk[0].length];
		for(int i=0;i<nnull;i++){
			for(int j=0;j<nxiz;j++) q_wk2[j][i] = 0;
		}
		for(int i=0;i<qq_wk[0].length;i++){
			for(int j=0;j<nxiz;j++) q_wk2[j][i+nnull] = qq_wk[j][i];
		}
		
		double[][] sr2 = new double[sr1.length+q_wk2.length][sr1[0].length];
		for(int i=0;i<sr1.length;i++) sr2[i] = sr1[i];
		for(int i=0;i<q_wk2.length;i++) sr2[i+sr1.length] = q_wk2[i];

		QRtest qr = new QRtest();
		QR sr = qr.qrtest(sr2);
		
		double[][] qrq = qr.qr_Q(sr);
		double trc = 0;
		for(int i=0;i<nobs;i++){
			for(int j=0;j<qrq[0].length;j++) trc += qrq[i][j]*qrq[i][j];
		}
		return trc;
	}
	
	public double[][] mkphi(Integer[] tm, double min, double max){
		double[][] s = new double[tm.length][2];
		for(int i=0;i<s.length;i++){
			s[i][0] = 1;
			s[i][1] = ((double)tm[i]-min)/(max-min)-0.5;
		}
		return s;
	}
	public double[][] mkrk(Integer[] tm, double[] x_basis, double min, double max){
		double[][] r = new double[tm.length][x_basis.length];
		double[] tm1 = new double[tm.length];
		for(int i=0;i<tm1.length;i++) tm1[i] = tm[i];
		for(int i=0;i<tm1.length;i++){
			tm1[i] = (tm1[i]-min)/(max-min);
		}
		for(int i=0;i<x_basis.length;i++){
			x_basis[i] = (x_basis[i]-min)/(max-min);
		}
		List<Double> Y = new ArrayList<Double>();
		List<Double> X = new ArrayList<Double>();
		for(int i=0;i<x_basis.length;i++){
			for(int j=0;j<tm1.length;j++){
				Y.add(x_basis[i]);
			}
		}
		int loop = (int) Math.ceil(Y.size()/tm1.length);
		for(int i=0;i<loop;i++){
			for(int j=0;j<tm1.length;j++){
				X.add(tm1[j]);
			}
		};
		double[] rk = new double[X.size()]; 
		for(int i=0;i<X.size();i++){
			double k2x = (Math.pow(X.get(i)-0.5,2)-(double)1/12)/2;
			double k2y = (Math.pow(Y.get(i)-0.5,2)-(double)1/12)/2;		
			double k4 = (Math.pow(Math.abs(X.get(i)-Y.get(i))-0.5,4) - Math.pow(Math.abs(X.get(i)-Y.get(i))-0.5,2)/2+(double)7/240)/24;
			rk[i] = k2x*k2y-k4;
		}
		int count = 0;
		for(int i=0;i<r[0].length;i++){
			for(int j=0;j<r.length;j++){
				r[j][i] = rk[count];
				count++;
			}
		}
		return r;
	}
	
	public double[][] transpose(double[][] a){
		double b[][] = new double[a[0].length][a.length];
		for(int i=1;i<=b.length;i++){
			for(int j=1;j<=b[0].length;j++){
				b[i-1][j-1] = a[j-1][i-1];
			}
		}
		return b;
	}
}
