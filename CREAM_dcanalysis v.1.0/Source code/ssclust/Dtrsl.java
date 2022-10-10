package ssclust;

public class Dtrsl {
	/*
	public static void main(String[] args) {
		double[] x = {1.1,2.1,3.1,4.1};

		int npar=4;
		Integer[] jpvt = {4,2,1,3};
		Integer job = 1;
		
		double[] result = dprmut(x, npar, jpvt, job);	
		System.out.println(result[0]+";"+result[1]+";"+result[2]+";"+result[3]);
	}
	*/
	public static Integer dtrsl(double[][] t, Integer ldt, Integer n, double[] b, Integer job, Integer info) {
		for(info=1;info<=n;info++){
			if(t[info-1][info-1] == 0){
				return info;
			}
		}
		info = 0;
		if(job == 01){
			b[n-1] = b[n-1]/t[n-1][n-1];
			if(n >= 2){
				for(int jj=2;jj<=n;jj++){
					int j = n-jj+1; //j从倒数第二个元素到第一个元素
					double temp = -b[j];
					double[] tt = new double[n];
					for(int m=1;m<=n;m++){
						tt[m-1] = t[m-1][j]; //tt为矩阵v的第j列
					}
					Reg.daxpy(j,temp,tt,1,b,1);
					b[j-1] = b[j-1]/t[j-1][j-1];
				}
			}
		}
		if(job == 11){
			b[0] = b[0]/t[0][0];
			if(n>=2){
				for(int j=2;j<=n;j++){ //j从第二个元素到最后一个元素
					double[] tt = new double[n]; 
					for(int m=1;m<=n;m++){
						tt[m-1] = t[m-1][j-1]; //tt为矩阵v的第j列
					}
					//System.out.println(j-1);
					//System.out.println(tt[0]+";"+tt[1]+";"+tt[2]+";"+tt[3]+";"+tt[4]+";"+tt[5]);
					//System.out.println(b[0]+";"+b[1]+";"+b[2]+";"+b[3]+";"+b[4]+";"+b[5]);
					//System.out.println("/////////////////");
					b[j-1] = b[j-1] - Reg.ddot(j-1,tt,b);
					b[j-1] = b[j-1]/t[j-1][j-1];
				}
			}
		}
		return info;
	}

}
