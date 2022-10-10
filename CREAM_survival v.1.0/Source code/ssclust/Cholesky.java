package ssclust;




public class Cholesky {
	/*
	public static void main(String[] args) {
		double[][] b = {{2.4,-0.8,1.4},{-0.8,-2.2,1.0},{1.4,1.0,2.3}};
		double[][] a = new double[b.length][b[0].length];
		for(int i=0;i<b.length;i++){
			for(int j=0;j<b[0].length;j++) a[i][j] = (double)b[i][j];
		}
		
		//double[][] a = {{2,-2},{0,5}};
		int lda=a.length;
		int p=a[0].length;
		double[] work = new double[p];
		Integer[] jpvt = new Integer[p];
		for(int i=0;i<jpvt.length;i++) jpvt[i] = 0;

		int info = a(a, lda, p, work, jpvt, 0);
		
		for(int i=0;i<a.length;i++) {
			for(int j=0;j<a[0].length;j++) {
				System.out.print(a[i][j]+"  ");
			}
			System.out.println();
		}
		System.out.println(info);
	}
	*/
	
	public static Integer a(double[][] a, Integer lda, Integer p, double[] work, Integer[] jpvt,
			Integer job) {
		/**cholesky·Ö½â**/
		int pl = 1;
		int pu = 0;
		int info = p;
		if(a.length != a[0].length){
			info = 0;
			return info;
		}
		for(int i=0;i<a[0].length;i++){
			for(int j=i+1;j<a.length;j++) a[j][i] = 0;
		}
		
		if(job !=0){
			for(int k=1;k<=p;k++){
				boolean swapk = jpvt[k-1]>0;
				boolean negk = jpvt[k-1]<0;
				jpvt[k-1] = k;
				if (negk) jpvt[k-1] = -jpvt[k-1];
				if(swapk){
					if(k!=pl){
						swap(pl-1,a,k-1,pl-1);
						double temp = a[k-1][k-1];
						a[k-1][k-1] = a[pl-1][pl-1];
						a[pl-1][pl-1] = temp;
						int plp1 = pl+1;
						if(p >= plp1){
							for(int j=plp1;j<=p;j++){
								if(j<k){
									double temp1 = a[pl-1][j-1];
									a[pl-1][j-1] = a[j-1][k-1];
									a[j-1][k-1] = temp1;
								}else{
									if(j>k){
									double temp1 = a[k-1][j-1];
									a[k-1][j-1] = a[pl-1][j-1];
									a[pl-1][j-1] = temp1;	
									}
								}
							}
						}
						jpvt[k-1] = jpvt[pl-1];
						jpvt[pl-1] = k;
					}
					pl=pl+1;
				}
			}
			pu = p;
			if(p>=pl){
				for(int kb=pl;kb<=p;kb++){
					int k=p-kb+pl;
					if(jpvt[k-1]<0){
						jpvt[k-1] = -jpvt[k-1];
						if(pu != k){
						swap(k-1,a,k-1,pu-1);
						double temp = a[k-1][k-1];
						a[k-1][k-1] = a[pu-1][pu-1];
						a[pu-1][pu-1] = temp;
						int kp1 = k + 1	;
						if(p>=kp1){
							for(int j=kp1;j<=p;j++){
								if(j<pu){
									double temp1 = a[k-1][j-1];
									a[k-1][j-1] = a[j-1][pu-1];
									a[j-1][pu-1] = temp1;
								}else{
									if(j>pu){
										double temp1 = a[k-1][j-1];
										a[k-1][j-1] = a[pu-1][j-1];
										a[pu-1][j-1] = temp1;
									}
								}
							}
						}
						int jt = jpvt[k-1];
						jpvt[k-1] = jpvt[pu-1];
						jpvt[pu-1] = jt;
						}
						pu--;
					}
				}
			}
		}
		//loop
		//for(int i=0;i<16;i++) jpvt[i] = 16-i;

		for(int k=1;k<=p;k++){		
			double maxdia = a[k-1][k-1];  		
			int kp1 = k + 1;
			int maxl = k;
			
			if(k>=pl && k<pu){
				for(int l=kp1;l<=pu;l++){
					if(a[l-1][l-1] > maxdia){
						maxdia = a[l-1][l-1];
						maxl = l;	
					}
				}
			}		
			if(maxdia<=0) {
				info = k-1;
				break;
			}
			
			if(k != maxl){
			int km1 = k - 1;
			swap(km1, a, k-1, maxl-1);
			a[maxl-1][maxl-1] = a[k-1][k-1];
			a[k-1][k-1] = maxdia;
			int jp = jpvt[maxl-1];
			jpvt[maxl-1] = jpvt[k-1];
			jpvt[k-1] = jp;
			}
			
			work[k-1] = Math.sqrt(a[k-1][k-1]);
			a[k-1][k-1] = work[k-1];
			if(p >= kp1){
				for(int j=kp1;j<=p;j++){
					if(k!=maxl){
						if(j<maxl){
							double temp = a[k-1][j-1];
							a[k-1][j-1] = a[j-1][maxl-1];
							a[j-1][maxl-1] = temp;
						}else{
							if(j>maxl){
								double temp = a[k-1][j-1];
								a[k-1][j-1] = a[maxl-1][j-1];
								a[maxl-1][j-1] = temp;
							}
						}
					}
					a[k-1][j-1] = a[k-1][j-1]/work[k-1];
					work[j-1] = a[k-1][j-1];
					double temp = -a[k-1][j-1];
					//a[kp1-1][j-1] = a[kp1-1][j-1] + temp*work[kp1-1];
					double[] dx = new double[j-k];
					for(int i=k+1;i<=j;i++){
						dx[i-k-1] = work[i-1];
					}
					double[] dy = new double[j-k];
					for(int i=k+1;i<=j;i++){
						dy[i-k-1] = a[i-1][j-1];
					}
					daxpy(j-k, temp, dx, 1, dy, 1);
					for(int i=k+1;i<=j;i++){
						a[i-1][j-1] = dy[i-k-1];
					}
				}
			}
		}
			
		return  info;
	}
	
	public static void swap(Integer p, double[][] a, Integer m, Integer n) {
		if(p>0){
			if(p<a.length && m<a[0].length && n<a[0].length){
				for(int i=0;i<p;i++){
					double t=a[i][m];
					a[i][m] = a[i][n];
					a[i][n] = t;
				}
			}else{
				System.out.println("Wrong size of p,m or n.");
			}
		}
		return;
	}
public static void daxpy (int n, double da, double dx[], int incx, double
            dy[], int incy) {
int i,ix,iy,m;

if (n <= 0) return;
if (da == 0.0) return;

if ((incx == 1) && (incy == 1)) {

//both increments equal to 1



for (i = 1; i <= n; i++) {
  dy[i-1] += da*dx[i-1];
}

return;
} else {
	
//at least one increment not equal to 1
ix = 1;
iy = 1;

if (incx < 0) ix = (-n+1)*incx + 1;
if (incy < 0) iy = (-n+1)*incy + 1;

for (i = 1; i <= n; i++) {

  dy[iy-1] += da*dx[ix-1];

  ix += incx;
  iy += incy;

}
return;
}
}
	
	
}
