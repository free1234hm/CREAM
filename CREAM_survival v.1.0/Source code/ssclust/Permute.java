package ssclust;

public class Permute {
	public static void main(String[] args) {
		double[] x = {1.1,2.1,3.1,4.1};

		int npar=4;
		Integer[] jpvt = {4,2,1,3};
		Integer job = 1;
		
		double[] result = dprmut(x, npar, jpvt, job);	
		System.out.println(result[0]+";"+result[1]+";"+result[2]+";"+result[3]);
	}
	public static double[] dprmut(double[] x, Integer npar, Integer[] jpvt, Integer job) {
		if(npar <= 1) return x;
		for(int j=1;j<=npar;j++){
			jpvt[j-1] = -jpvt[j-1];
		}
		if(job == 0){
			for(int i=1;i<=npar;i++){
				if(jpvt[i-1] <=0){
					int j=i;
					jpvt[j-1] = -jpvt[j-1];
					int k = jpvt[j-1];
					while(jpvt[k-1] < 0){
						double t = x[j-1];
						x[j-1] = x[k-1];
						x[k-1] = t;
						jpvt[k-1] = -jpvt[k-1];
						j=k;
						k=jpvt[k-1];
					}
				}
			}
		}else{
			for(int i=1;i<=npar;i++){
				if(jpvt[i-1] <=0){
					jpvt[i-1] = -jpvt[i-1];
				    int j = jpvt[i-1];
				    while(j != i){
				    	double t = x[i-1];
						x[i-1] = x[j-1];
						x[j-1] = t;
						jpvt[j-1] = -jpvt[j-1];
					    j = jpvt[j-1];
				    }
				}
			}
		}
		
		return x;
	}
	
}
