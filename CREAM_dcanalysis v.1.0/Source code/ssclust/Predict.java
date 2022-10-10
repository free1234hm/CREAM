package ssclust;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.List;

public class Predict {
	public List<Object> predict(Ssan object, List<Object> newdata, boolean se_fit) {		
		List<Object> result = new ArrayList<Object>();
		int nnull, nz;
		int nnew = ((double[])newdata.get(0)).length;
	    int nbasis = object.id_basis.length;
	    if(object.z.d != null && object.z.d.length > 0) {
	    	nnull = object.z.d.length;
	    }else{
	    	nnull = 0;
	    }
	    if(object.z.b != null && object.z.b.length > 0) {
	    	nz = object.z.b.length;
	    	
	    }else{
	    	nz = 0;
	    }
	    
	    List<Object> term = object.term;
	    double env_min = ((Labels)term.get(1)).env_min;
	    double env_max = ((Labels)term.get(1)).env_max;
	    
	    int nq = 0;
	    double[][] s = new double[nnew][term.size()];
	    double[][] r = new double[nnew][nbasis];
	    List<Integer> philist = new ArrayList<Integer>();
	    List<Integer> rklist = new ArrayList<Integer>();
	    for(int label=1;label<=term.size();label++){
	    	if (label==1) {
	    		philist.add(((Labels)term.get(label-1)).iphi);
	    		for(int i=1;i<=nnew;i++) s[i-1][label-1] = 1;
	    		continue;
	    	}
	    	if(label == 3) continue;
	    	double[] xnew = (double[])newdata.get(0);
	    	double[] x = new double[object.id_basis.length]; 	
	    	for(int i=1;i<=x.length;i++) x[i-1] = ((Integer[])object.mf.get(label-1))[object.id_basis[i-1]-1];
	    	int nphi = ((Labels)term.get(label-1)).nphi;
	        int nrk = ((Labels)term.get(label-1)).nrk;
	        if(nphi>0){
	        	int iphi = ((Labels)term.get(label-1)).iphi;
	            for(int i=1;i<=nphi;i++){
	            	philist.add(iphi+(i-1));
	            	for(int j=1;j<=s.length;j++){
	        			s[j-1][label-1] = (xnew[j-1]-env_min)/(env_max-env_min)-0.5;
	        		}
	            }
	        }
	        
	        if(nrk>0){
	        	int irk = ((Labels)term.get(label-1)).irk;
	        	 for(int i=1;i<=nrk;i++){
	        		 rklist.add(irk+(i-1)); 
	        		 nq=nq+1;
	        		 //System.out.println(xnew[0]+"  "+xnew[1]+"  "+xnew[2]+"  "+xnew[3]+"  "+xnew[4]);
	        		 //System.out.println(x[0]+"  "+x[1]+"  "+x[2]+"  "+x[3]+"  "+x[4]);
	        		 r = mkrk(xnew, x, env_min, env_max);
	        	 }
	        } 
	    }
	    
	    double[][] r_wk1 = new double[nnew][nbasis]; 
	    for(int i=1;i<=rklist.size();i++){
	    	for(int m=1;m<=r_wk1.length;m++){
	    		for(int n=1;n<=r_wk1[0].length;n++){
	    			r_wk1[m-1][n-1] = Math.pow(10, object.z.theta)*r[m-1][n-1];
	    		}
	    	}
	    }
	    
	    double[][] r_wk = null; 
	    if(nz>0){
	    	r_wk = new double[nnew][nbasis+nz];
	    	for(int m=1;m<=r_wk1.length;m++){
	    		for(int n=1;n<=r_wk1[0].length;n++){
	    			r_wk[m-1][n-1] = r_wk1[m-1][n-1];
	    		}
	    	}
	    }else{
	    	r_wk = r_wk1;
	    }
	    
	    
	    
	    
	    int nphi = philist.size();
	    double[][] cb = null;
	    if(object.z.b != null && object.z.b.length > 0){
	    	cb = new double[object.z.c.length+object.z.b.length][1];
		    for(int i=1;i<=object.z.c.length;i++) cb[i-1][0] = object.z.c[i-1];
		    for(int i=1;i<=object.z.b.length;i++) cb[i+object.z.c.length-1][0] = object.z.b[i-1];
		    
	    }else{
	    	cb = new double[object.z.c.length][1];
		    for(int i=1;i<=cb.length;i++) cb[i-1][0] = object.z.c[i-1];
	    }
	    
	   
	    double[][] pmean1 = Matrixdot(r_wk,cb); 
	    /*
	    for(int i=0;i<pmean1.length;i++){
	    	for(int j=0;j<pmean1[0].length;j++){
	    		System.out.print(pmean1[i][j]+" ");
	    	}
	    	System.out.println();
	    }
	    System.out.println("////////////////////////");
	    */
	    if (nphi>0){ 
	    	double[][] dd = new double[philist.size()][1];
	    	for(int i=1;i<=dd.length;i++) dd[i-1][0] = object.z.d[philist.get(i-1)-1];
	    	double[][] pmean2 = Matrixdot(s,dd);
	    	for(int i=1;i<=pmean1.length;i++){
	    		for(int j=1;j<=pmean1[0].length;j++){
	    			pmean1[i-1][j-1] = pmean1[i-1][j-1] + pmean2[i-1][j-1];
	    		}
	    	}
	    }
	    
	    /*
	    double[][] wk1 = new double[r_wk.length][r_wk[0].length+s[0].length];
	    for(int i=0;i<wk1.length;i++){
	    	for(int j=0;j<s[0].length;j++){
	    		wk1[i][j] = s[i][j];
	    	}
	    	for(int j=0;j<r_wk[0].length;j++){
	    		wk1[i][j+s[0].length] = r_wk[i][j];
	    	}
	    }
	    double[][] dd = new double[philist.size()][1];
    	for(int i=1;i<=dd.length;i++) dd[i-1][0] = object.z.d[philist.get(i-1)-1];
    	double[][] cbd = new double[cb.length+dd.length][cb[0].length];
    	 for(int i=0;i<dd.length;i++){
    		 cbd[i][0] = dd[i][0];
    	 }
    	 for(int i=0;i<cb.length;i++){
    		 cbd[i+dd.length][0] = cb[i][0];
    	 }
    	 double[][] mean = Matrixdot(wk1,cbd);
	    */
	    
	    
	    double[] pmean = new double[pmean1.length];
	    for(int i=1;i<=pmean1.length;i++) pmean[i-1] = pmean1[i-1][0];
	    
	    if(object.term.size()>2){
	    	double[] offset;
	    	if(newdata.size() == 1){
	    		javax.swing.JOptionPane.showMessageDialog(null, "missing offset!");
				System.exit(0);
	    	}
	    	offset = (double[])newdata.get(1);
	    	for(int i=0;i<pmean.length;i++) pmean[i] = pmean[i] + offset[i];
	    }
	    
	    result.add(pmean);
	    
	    if(se_fit == true){
	    	double b = object.z.fit.varht/Math.pow(10, object.z.nlambda);
	    	double[][] ss = new double[nnull][nnew];
	    	if (philist != null && philist.size()>0){
	    		double[][] ts = transpose(s);
	    		for(int i=1;i<=philist.size();i++) ss[philist.get(i-1)-1] = ts[i-1];
	    	    }
	    		double[][] r2 = Matrixdot(r_wk, object.z.se_aux_ssanova.vec);
	    		double[][] rr = transpose(r2);
	    		double[][] sr = new double[ss.length+rr.length][ss[0].length];
	    		for(int i=1;i<=ss.length;i++) sr[i-1] = ss[i-1];
	    		for(int i=1;i<=ss.length;i++) sr[i+ss.length-1] = rr[i-1];
	    		double[][] wk = Matrixdot(object.z.se_aux_ssanova.hfac, sr);
	    		double[] pse = new double[wk[0].length];
	    		for(int i=1;i<=wk[0].length;i++){
	    			for(int j=1;j<=wk.length;j++) pse[i-1] = pse[i-1] + wk[j-1][i-1]*wk[j-1][i-1];
	    			pse[i-1] = b*pse[i-1];
	    			pse[i-1] = Math.sqrt(pse[i-1]);
	    		}
	    		result.add(pse);
	    }    	
	    return result;   
	}
	
	
	
	static double[][] Matrixdot(double[][] a,double[][] b){
		double c[][] = new double[a.length][b[0].length];
		int x,i,j;
		for(i = 0;i<a.length ;i++)
		{
			for(j = 0;j<b[0].length;j++)
			{
				double temp = 0;
				for(x = 0;x<b.length;x++)
				{
					temp+=a[i][x]*b[x][j];
					
				}
				c[i][j] = temp;
				
			}
		}
		return c;
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
	public double[][] mkphi(double[] tm, double min, double max){
		double[][] s = new double[tm.length][2];
		for(int i=0;i<s.length;i++){
			s[i][0] = 1;
			s[i][1] = (tm[i]-min)/(max-min)-0.5;
		}
		return s;
	}
	public double[][] mkrk(double[] tm, double[] x_basis, double min, double max){
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

}
