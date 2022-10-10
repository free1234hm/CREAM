package ssclust;

public class QRtest {

	public static void main(String[] args) {
		
		double[][] aa = {{1,2,5},{2,3,8},{6,7,21}};
		double[][] y = {{1,2},{2,3},{6,7}};
		QRtest qr = new QRtest();
		QR wk1 = qr.qrtest(aa);
		//System.out.println(wk1.qr[0][0]+";"+wk1.qr[0][1]+";"+wk1.qr[0][2]);
		double[][] qy = qr.qrrsd(wk1, y);
		//System.out.println(wk1.qr[0][0]+";"+wk1.qr[0][1]+";"+wk1.qr[0][2]);
		//System.out.println(wk1.qr[1][0]+";"+wk1.qr[1][1]+";"+wk1.qr[1][2]);
		//System.out.println(wk1.qr[2][0]+";"+wk1.qr[2][1]+";"+wk1.qr[2][2]);
		//System.out.println(wk1.qraux[0]+";"+wk1.qraux[1]+";"+wk1.qraux[2]);
		//System.out.println(wk1.rank);
		//System.out.println(wk1.pivot[0]+";"+wk1.pivot[1]+";"+wk1.pivot[2]);
		System.out.println(qy[0][0]+";"+qy[0][1]);
		System.out.println(qy[1][0]+";"+qy[1][1]);
		System.out.println(qy[2][0]+";"+qy[2][1]);
	}
	
	
	public QR qrtest(double[][] x){
		  QR result = new QR(); 
		  int p = x[0].length;
		  int n = x.length;  
		  int ldx = n;
		  double tol = 1e-07;
		  int rank;
		  Integer[] jpvt = new Integer[p];
		  for(int i=0;i<jpvt.length;i++) jpvt[i] = i+1;
		  double[] qraux = new double[p];
		  double[][] work = new double[p][2];
		  
		  for(int j=1;j<=p;j++){
			  double[] xx = new double[n];
			  for(int i=1;i<=n;i++){
				  xx[i-1] = x[i-1][j-1];
			  }
			  qraux[j-1] = dnrm2(n,xx,1);
			  work[j-1][0] = qraux[j-1];
			  work[j-1][1] = qraux[j-1];
			  if(work[j-1][1] == 0) work[j-1][1] = 1;
		  }
		  int lup = Math.min(n,p);
		  int k = p + 1;
		  
		  for(int l=1;l<=lup;l++){
			  while(l<k && qraux[l-1] < work[l-1][1]*tol){
				  int lp1 = l+1;
				  for(int i=1;i<=n;i++){
					  double t = x[i-1][l-1];
					  for(int j=lp1;j<=p;j++){
						  x[i-1][j-2] = x[i-1][j-1];
					  }
					  x[i-1][p-1] = t;
				  }
				  int i = jpvt[l-1];
				  double t = qraux[l-1];
				  double tt = work[l-1][0];
				  double ttt = work[l-1][1];
		          for(int j=lp1;j<=p;j++){
		        	  jpvt[j-2] = jpvt[j-1];
		        	  qraux[j-2] = qraux[j-1];
		        	  work[j-2][0] = work[j-1][0];
		        	  work[j-2][1] = work[j-1][1];
		          }
		          jpvt[p-1] = i;
		          qraux[p-1] = t;
		          work[p-1][0] = tt;
		          work[p-1][1] = ttt;
		          k = k - 1;
			  }
			  if(l != n){
				  double[] xx = new double[n-l+1];
				  for(int i=l;i<=n;i++){
					  xx[i-l] = x[i-1][l-1];
				  }
				  double nrmxl = dnrm2(n-l+1,xx,1);
				  if(nrmxl != 0){
					  if (x[l-1][l-1] != 0){
						  if(x[l-1][l-1]<0) nrmxl = -nrmxl;
					  }
					  for(int i=l;i<=n;i++){
						  x[i-1][l-1] = x[i-1][l-1]/nrmxl;
					  }
					  x[l-1][l-1] = 1 + x[l-1][l-1];
					  int lp1 = l + 1;
					  if(p>=lp1){
						  for(int j=lp1;j<=p;j++){
							  double[] x1 = new double[n-l+1];
							  for(int i=l;i<=n;i++){
								  x1[i-l] = x[i-1][l-1];
							  }
							  double[] x2 = new double[n-l+1];
							  for(int i=l;i<=n;i++){
								  x2[i-l] = x[i-1][j-1];
							  }
							  double t = -Reg.ddot(n-l+1,x1,x2)/x[l-1][l-1];
							  Reg.daxpy(n-l+1,t,x1,1,x2,1);
							  
							  for(int i=l;i<=n;i++){
								   x[i-1][j-1] = x2[i-l];
							  }
							  
							  if(qraux[j-1] !=0){
								  double tt = 1 - (Math.abs(x[l-1][j-1])/qraux[j-1])*(Math.abs(x[l-1][j-1])/qraux[j-1]);
						          tt = Math.max(tt,0);
						          t = tt;
						          if(Math.abs(t) >= 1e-06){
						        	  qraux[j-1] = qraux[j-1]*Math.sqrt(t);
						          }else{
						        	  double[] x3 = new double[n-l];
									  for(int i=l+1;i<=n;i++){
										  x3[i-l-1] = x[i-1][l-1];
									  }
						        	  qraux[j-1] = dnrm2(n-l,x3,1);
						              work[j-1][0] = qraux[j-1];
						          }
							  }
						  }
					  }
					  qraux[l-1] = x[l-1][l-1];
					  x[l-1][l-1] = -nrmxl;
				  }
			  }
			  
		  }
		  k = Math.min(k - 1, n); 
		  
		  result.qr = x;
		  result.qraux = qraux;
		  result.pivot = jpvt;
		  result.rank = k;
		  return result;
	}
	
	public double[][] qrqty(QR qr, double[][] y){
		int n = qr.qr.length;
	    int k = qr.rank;
	    int ny = y[0].length;
	    if(y.length != n) {
	    	javax.swing.JOptionPane.showMessageDialog(null, "'qr' and 'y' must have the same number of rows");
	    	System.exit(0);
	    }
	    double[][] qty = y;
	    dqrqty(qr.qr, n, k, qr.qraux, y, ny, qty);
	    return qty; 
	}
	
	public void dqrqty(double[][] x, Integer n, Integer k, double[] qraux, double[][] y, Integer ny,
			double[][] qty){
		int info = 0;
		for(int j=1;j<=ny;j++){
			double[] yy = new double[y.length];
			for(int i=1;i<=yy.length;i++) yy[i-1] = y[i-1][j-1];
			double[] qtyy = new double[y.length];
			for(int i=1;i<=qtyy.length;i++) qtyy[i-1] = qty[i-1][j-1];
			double[] qy = qtyy;
			double[] rsd = qtyy;
			dqrsl(x, n, n, k, qraux, yy, qtyy, qy, rsd, info, "cqty");
			for(int i=1;i<=qtyy.length;i++) qty[i-1][j-1] = qtyy[i-1];
		}
	}
	
	public double[][] qr_Q(QR qr){
		Integer[] dqr = {qr.qr.length, qr.qr[0].length};
		int n = dqr[0];
		int min;
		if(dqr[0]<dqr[1]){
			min = dqr[0];
		}else{
			min = dqr[1];
		}
		double[] Dvec = new double[min];
		for(int i=0;i<Dvec.length;i++) Dvec[i] = 1;
		double[][] D = new double[n][min];
		for(int i=0;i<min;i++) D[i][i] = 1;
		double[][] qy = dqrqy(qr, D);
		return qy;
	}
	
	public double[][] dqrqy(QR qr, double[][] y){
		int n = qr.qr.length;
		int k = qr.rank;
		int ny = y[0].length;
		int info = 0;
		double[][] qy = y;
		for(int j=1;j<=ny;j++){
			double[] yy = new double[y.length];
			for(int i=1;i<=yy.length;i++) yy[i-1] = y[i-1][j-1];
			double[] qyy = new double[y.length];
			for(int i=1;i<=qyy.length;i++) qyy[i-1] = qy[i-1][j-1];
			double[] qty = qyy;
			double[] rsd = qyy;
			dqrsl(qr.qr, n, n, k, qr.qraux, yy, qty, qyy, rsd, info, "cqy");
			for(int i=1;i<=qyy.length;i++) qy[i-1][j-1] = qyy[i-1];
		}
		return qy;
	}
	
	public double[][] qrrsd(QR qr, double[][] y){
		int n = qr.qr.length;
		//System.out.println(n);
	    int k = qr.rank;
	    int ny = y[0].length;
	    if(y.length != n) {
	    	javax.swing.JOptionPane.showMessageDialog(null, "'qr' and 'y' must have the same number of rows");
	    	System.exit(0);
	    }
	    double[][] rsd = y;
	    dqrrsd(qr.qr, n, k, qr.qraux, y, ny, rsd);
	    return rsd; 
	}
	
	public void dqrrsd(double[][] x, Integer n, Integer k, double[] qraux, double[][] y, Integer ny,
			double[][] rsd){
		int info = 0;
		for(int j=1;j<=ny;j++){
			double[] yy = new double[y.length];
			for(int i=1;i<=yy.length;i++) yy[i-1] = y[i-1][j-1];
			double[] rsdd = new double[y.length];
			for(int i=1;i<=rsdd.length;i++) rsdd[i-1] = rsd[i-1][j-1];
			double[] qty = rsdd;
			double[] qy = rsdd;
			dqrsl(x, n, n, k, qraux, yy, qty, qy, rsdd, info, "cr");
			for(int i=1;i<=rsdd.length;i++) rsd[i-1][j-1] = rsdd[i-1];
		}
	}
	
	public void dqrsl(double[][] x, Integer ldx, Integer n, Integer k, double[] qraux, double[] y,
			 double[] qty, double[] qy, double[] rsd, Integer info, String method){
		info=0;
		int ju = Math.min(k,n-1);
		if(ju == 0){
			qty[0] = y[0];
			return;
		}
		if(method.equals("cqty") || method.equals("cr")){
			Reg.dcopy(n,y,1,qty,1);
			for(int j=1;j<=ju;j++){
				if(qraux[j-1] != 0){
					double temp = x[j-1][j-1];
					x[j-1][j-1] = qraux[j-1];
					double[] x1 = new double[n-j+1];
					for(int m=j;m<=n;m++){
						x1[m-j] = x[m-1][j-1];
					}
					double[] x2 = new double[n-j+1];
					for(int m=j;m<=n;m++){
						x2[m-j] = qty[m-1];
					}
					double t = -Reg.ddot(n-j+1,x1,x2)/x[j-1][j-1];
					Reg.daxpy(n-j+1,t,x1,1,x2,1);
					for(int m=j;m<=n;m++){
						qty[m-1] = x2[m-j];
					}
					x[j-1][j-1] = temp;
				}
			}
		}
		if(method.equals("cqy")){
			Reg.dcopy(n,y,1,qy,1);
			for(int jj=1;jj<=ju;jj++){
				int j = ju - jj + 1;
				if(qraux[j-1] != 0){
					double temp = x[j-1][j-1];
					x[j-1][j-1] = qraux[j-1];
					double[] x1 = new double[n-j+1];
					for(int m=j;m<=n;m++){
						x1[m-j] = x[m-1][j-1];
					}
					double[] x2 = new double[n-j+1];
					for(int m=j;m<=n;m++){
						x2[m-j] = qy[m-1];
					}
					double t = -Reg.ddot(n-j+1,x1,x2)/x[j-1][j-1];
					Reg.daxpy(n-j+1,t,x1,1,x2,1);
					for(int m=j;m<=n;m++){
						qy[m-1] = x2[m-j];
					}
					x[j-1][j-1] = temp;
				}
			}
		}
		int kp1 = k+1;
		if(method.equals("cr") && k<n){	
			double[] qty2 = new double[n-k];
			for(int i=kp1;i<=n;i++) qty2[i-k-1] = qty[i-1];
			double[] rsd2 = new double[n-k];
			for(int i=kp1;i<=n;i++) rsd2[i-k-1] = rsd[i-1];
			Reg.dcopy(n-k,qty2,1,rsd2,1);
			for(int i=kp1;i<=n;i++) rsd[i-1] = rsd2[i-k-1];
		}
		if(method.equals("cr")){
			for(int i=1;i<=k;i++) rsd[i-1] = 0;
		}
		if(method.equals("cr")){
			for(int jj=1;jj<=ju;jj++){
				int j = ju - jj + 1;
				if(qraux[j-1] != 0){
					double temp = x[j-1][j-1];
					x[j-1][j-1] = qraux[j-1];
					double[] x1 = new double[n-j+1];
					for(int m=j;m<=n;m++){
						x1[m-j] = x[m-1][j-1];
					}
					double[] x2 = new double[n-j+1];
					for(int m=j;m<=n;m++){
						x2[m-j] = rsd[m-1];
					}
					double t = -Reg.ddot(n-j+1,x1,x2)/x[j-1][j-1];
					Reg.daxpy(n-j+1,t,x1,1,x2,1);
					for(int m=j;m<=n;m++){
						rsd[m-1] = x2[m-j];
					}
					x[j-1][j-1] = temp;
				}
			}
		}
		
	}
	
	public double dnrm2(int n, double[] x, int incx){
		double sum = 0;
		double result = 0;
		for(int i=1;i<=n;i++){
			sum = sum + x[i-1]*x[i-1];
		}
		result = Math.sqrt(sum);
		return result;
	}
	
}
