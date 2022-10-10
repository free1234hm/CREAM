package ssclust;

public class Reg {
	public static Regression reg(int run, double[][] sr, Integer nobs, Integer nnull, double[][] q, Integer nxi,
			double[] y, double alpha, double varht, double score, double[] dc, double mchpr,
			double[][] v, double[] mu, Integer[] jpvt, double[] wk, Integer rkv) {
		
		Regression result = new Regression();
		int nn = nnull + nxi; //nnull=2,nxi=时间点数+类基因数
		
		int infowk = 0;
		for (int i=1;i<=nn;i=i+1){
			infowk = infowk + jpvt[i-1];
		}
		
		rkv = Cholesky.a(v, nn, nn, wk, jpvt, 1); //Cholesky decomposition
		double[] idamaxv = new double[nn*nn];
		int count = 0;
		for(int i=infowk+1;i<=nn;i++){
			idamaxv[count] = v[i-1][infowk];
			count++;
		}
		for(int i=infowk+2;i<=nn;i++){
			for(int j=1;j<=nn;j++){
			idamaxv[count] = v[j-1][i-1];
			count++;
			}
		}
		int jj = idamax(rkv-infowk, idamaxv, nn+1);
		while(v[rkv-1][rkv-1]<v[infowk+jj-1][infowk+jj-1]*Math.sqrt(mchpr)){
			rkv--;
		}
		
		for(int i=rkv+1;i<=nn;i++){
			if(jj>0){
				v[i-1][i-1] = v[jj-1][jj-1];	
			}else{
				v[i-1][i-1] = 1.240105e-321;
			}			
			
			double[] idsetv = new double[nn*nn];
			int count2 = 0;
			for(int j=rkv+1;j<=nn;j++){
				idsetv[count2] = v[j-1][i-1];
				count2++;
			}
			for(int j=i+1;j<=nn;j++){
				for(int m=1;m<=nn;m++){
				idsetv[count2] = v[m-1][j-1];
				count2++;
				}
			}
			dset(i-rkv-1, 0.0, idsetv, 1);
			int count3 = 0;
			for(int j=rkv+1;j<=nn;j++){
				v[j-1][i-1] = idsetv[count3];
				count3++;
			}
			for(int j=i+1;j<=nn;j++){
				for(int m=1;m<=nn;m++){
				v[m-1][j-1] = idsetv[count3];
				count3++;
				}
			}
		}
		
		dcopy (nn, mu, 1, dc, 1); //mu已知，所以此步dc已知
		dc = Permute.dprmut (dc, nn, jpvt, 0);		
		
		infowk = Dtrsl.dtrsl (v, nn, nn, dc, 11, infowk); //dc变化，v不变
		
		
		double[] idsetv = new double[nn];
		for(int i=rkv+1;i<=nn;i++){
			idsetv[i-rkv-1] = dc[i-1];
		}
		dset (nn-rkv, 0.0, idsetv, 1);
		for(int i=rkv+1;i<=nn;i++){
			dc[i-1] = idsetv[i-rkv-1]; 
		}
		infowk = Dtrsl.dtrsl (v, nn, nn, dc, 01, infowk);
		
		dc = Permute.dprmut (dc, nn, jpvt, 1);
		
		for (int i=1;i<=nobs;i++){
			double[] sr_row = new double[nn];
			for(int j=1;j<=nn;j++){
				sr_row[j-1]=sr[i-1][j-1];
			}
			wk[i-1] = y[i-1] - ddot(nn, sr_row, dc); //各基因在各时间点表达水平-各时间点中心
		}
		double rss=0; 
		double trc=0;

			rss = ddot (nobs, wk, wk) / nobs;	
			for (int i=1;i<=nobs;i++) {
				double[] sr_row = new double[nn];
				for(int j=1;j<=nn;j++){
					sr_row[j-1]=sr[i-1][j-1];
				}
		        dcopy (nn, sr_row, 1, mu, 1);
		        Permute.dprmut (mu, nn, jpvt, 0);   
		        Dtrsl.dtrsl (v, nn, nn, mu, 11, infowk);
		        wk[i-1] = ddot (nn, mu, mu);
		    }

			double dasum=0;
			for(int i=1;i<=nobs;i++){
				dasum = dasum + Math.abs(wk[i-1]);
			}
			trc = dasum / nobs;

		    score = rss / ((1.0-alpha*trc)*(1.0-alpha*trc));
		    varht = rss / (1.0-trc);

		wk[0] = rss;
		wk[1] = trc;
		
		result.score = score;
		result.varht = varht;
		result.dc = dc;
		result.chol = v;
		result.jpvt = jpvt;
		result.wk = wk;
		result.rkv = rkv;
		return  result;
	}
	
	public static Regression reg(int run, double[][] sr, Integer nobs, Integer nnull, double[][] q, Integer nxi,
			double[] y, Integer method, double alpha, double varht, double score, double[] dc, double mchpr,
			double[][] v, double[] mu, Integer[] jpvt, double[] wk, Integer rkv) {
		
		Regression result = new Regression();
		int nn = nnull + nxi; //nnull=2,nxi=时间点数+类基因数
		
		int infowk = 0;
		for (int i=1;i<=nn;i=i+1){
			infowk = infowk + jpvt[i-1];
		}
		
		rkv = Cholesky.a(v, nn, nn, wk, jpvt, 1); //Cholesky decomposition
		double[] idamaxv = new double[nn*nn];
		int count = 0;
		for(int i=infowk+1;i<=nn;i++){
			idamaxv[count] = v[i-1][infowk];
			count++;
		}
		for(int i=infowk+2;i<=nn;i++){
			for(int j=1;j<=nn;j++){
			idamaxv[count] = v[j-1][i-1];
			count++;
			}
		}
		int jj = idamax(rkv-infowk, idamaxv, nn+1);
		while(v[rkv-1][rkv-1]<v[infowk+jj-1][infowk+jj-1]*Math.sqrt(mchpr)){
			rkv--;
		}
		
		for(int i=rkv+1;i<=nn;i++){
			if(jj>0){
				v[i-1][i-1] = v[jj-1][jj-1];	
			}else{
				v[i-1][i-1] = 1.240105e-321;
			}			
			
			double[] idsetv = new double[nn*nn];
			int count2 = 0;
			for(int j=rkv+1;j<=nn;j++){
				idsetv[count2] = v[j-1][i-1];
				count2++;
			}
			for(int j=i+1;j<=nn;j++){
				for(int m=1;m<=nn;m++){
				idsetv[count2] = v[m-1][j-1];
				count2++;
				}
			}
			dset(i-rkv-1, 0.0, idsetv, 1);
			int count3 = 0;
			for(int j=rkv+1;j<=nn;j++){
				v[j-1][i-1] = idsetv[count3];
				count3++;
			}
			for(int j=i+1;j<=nn;j++){
				for(int m=1;m<=nn;m++){
				v[m-1][j-1] = idsetv[count3];
				count3++;
				}
			}
		}
		
			

		dcopy (nn, mu, 1, dc, 1); //mu已知，所以此步dc已知
		dc = Permute.dprmut (dc, nn, jpvt, 0);
		
		
		infowk = Dtrsl.dtrsl (v, nn, nn, dc, 11, infowk); //dc变化，v不变
		
		
		double[] idsetv = new double[nn];
		for(int i=rkv+1;i<=nn;i++){
			idsetv[i-rkv-1] = dc[i-1];
		}
		dset (nn-rkv, 0.0, idsetv, 1);
		for(int i=rkv+1;i<=nn;i++){
			dc[i-1] = idsetv[i-rkv-1]; 
		}
		infowk = Dtrsl.dtrsl (v, nn, nn, dc, 01, infowk);
		
		dc = Permute.dprmut (dc, nn, jpvt, 1);
		
		
		
			
		
		for (int i=1;i<=nobs;i++){
			double[] sr_row = new double[nn];
			for(int j=1;j<=nn;j++){
				sr_row[j-1]=sr[i-1][j-1];
			}
			wk[i-1] = y[i-1] - ddot(nn, sr_row, dc); //各基因在各时间点表达水平-各时间点中心
		}
		double rss=0; 
		double trc=0;

			rss = ddot (nobs, wk, wk) / nobs;	
			for (int i=1;i<=nobs;i++) {
				double[] sr_row = new double[nn];
				for(int j=1;j<=nn;j++){
					sr_row[j-1]=sr[i-1][j-1];
				}
		        dcopy (nn, sr_row, 1, mu, 1);
		        Permute.dprmut (mu, nn, jpvt, 0);   
		        Dtrsl.dtrsl (v, nn, nn, mu, 11, infowk);
		        wk[i-1] = ddot (nn, mu, mu);
		    }

			double dasum=0;
			for(int i=1;i<=nobs;i++){
				dasum = dasum + Math.abs(wk[i-1]);
			}
			trc = dasum / nobs;
			
			if (method==2) {
		        score = rss / ((1.0-alpha*trc)*(1.0-alpha*trc));
		        varht = rss / (1.0-trc);
		    }else{
		    	score = rss + 2.0 * varht * alpha * trc;
		    }

		wk[0] = rss;
		wk[1] = trc;
		
		result.score = score;
		result.varht = varht;
		result.dc = dc;
		result.chol = v;
		result.jpvt = jpvt;
		result.wk = wk;
		result.rkv = rkv;
		return  result;
	}
	
	public static int idamax (int n, double x[], int incx) {
	      double xmax;
	      int isamax,i,ix;
	      if (n < 1 || incx <=0) {
	         isamax = 0;
	         return isamax;
	      } else if (n == 1) {
	         isamax = 1;
	         return isamax;
	      } else if (incx == 1) {
	         isamax = 1;
	         xmax = Math.abs(x[0]);         
	         for (i = 2; i <= n; i++) {
	            if (Math.abs(x[i-1]) > xmax) {
	               isamax = i;
	               xmax = Math.abs(x[i-1]);
	            }
	         }
	      } else {
	         isamax = 1;
	         ix = 1;
	         xmax = Math.abs(x[ix-1]);
	         ix += incx;
	         for (i = 2; i <= n; i++) {
	            if (Math.abs(x[ix-1]) > xmax) {
	               isamax = i;
	               xmax = Math.abs(x[ix-1]);
	            }
	            ix += incx;
	         }	        
	      }
	      return isamax;
	   }
	
	public static double ddot(Integer n, double[] x, double[] y) {
		double ddot = 0;
		if(n<=0){
			return ddot;
		}
		for(int i=1;i<=n;i++){
			ddot=ddot+x[i-1]*y[i-1];
		}
		return ddot;
	}
	
	public static void dcopy(Integer n, double[] x, Integer incx, double[] y, Integer incy) {
		if(n<=0){
			return;
		}
		if(incx ==1 && incy == 1){
			for(int i=1;i<=n;i++){
				y[i-1] = x[i-1];
			}
			return;
		}
		
	}
	
	public static void dset(Integer n, double alpha, double[] x, Integer incx) {
		if(n<=0) return;
		if(incx != 1){
			int ix;
			if(incx>0){
				ix=1;
			}else{
				ix=1-(n-1)*incx;
			}
			for(int i=1;i<=n;i++){
				x[ix-1]=alpha;
				ix=ix+incx;
			}
			return;
		}else{
			for(int i=1;i<=n;i++){
				x[i-1] = alpha;
			}
			return;
		}
	}
	
	public static void daxpy (int n, double da, double dx[], int incx, double
            dy[], int incy) {
           int i,ix,iy,m;
           if (n <= 0) return;
           if (da == 0.0) return;
           if ((incx == 1) && (incy == 1)) {
               //both increments equal to 1
               m = n%4;
               for (i = 1; i <= m; i++) {
                   dy[i-1] += da*dx[i-1];
                   }
               for (i = m+1; i <= n; i += 4) {
                   dy[i-1]   += da*dx[i-1];
                   dy[i] += da*dx[i];
                   dy[i+1] += da*dx[i+1];
                   dy[i+2] += da*dx[i+2];
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
