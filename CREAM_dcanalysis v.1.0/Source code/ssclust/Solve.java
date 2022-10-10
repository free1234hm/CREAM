package ssclust;

public class Solve {
	public double[] caculate(double[][] a)// 高斯消元求未知数X，     
    {          
    int _rows = a.length;         
    int _cols = a[0].length;         
    int L = _rows - 1;         
    int i, j, l, n, m, k = 0;          
    double[] temp1 = new double[_rows];          /* 第一个do-while是将增广矩阵消成上三角形式 */         
    do {    
              n = 0;
              for (l = k; l < L; l++)  temp1[n++] = a[l + 1][k] / a[k][k];             
              for (m = 0, i = k + 1; i < _rows; i++, m++) {                 
              for (j = k; j < _cols; j++)  a[i][j] -= temp1[m] * a[k][j]; 
               }            
                k++;        
                 } while (k < _rows);          // /*第二个do-while是将矩阵消成对角形式，并且重新给k赋值,最后只剩下对角线和最后一列的数，其它都为0*/         
                 k = L - 1;         
                 do {             
                  n = 0;              
                  for (l = k; l >= 0; l--) temp1[n++] = a[k - l][k + 1] / a[k + 1][k + 1];             
                  for (m = 0, i = k; i >= 0; i--, m++) {                 
                  for (j = k; j < _cols; j++)  a[k - i][j] -= temp1[m] * a[k + 1][j];
                   }             
                   k--;         
                   }while (k >= 0); 
       double[] newresult = new double[_rows];         
       for (i = 0; i < _rows; i++) {              
       newresult[i] = a[i][_rows] / a[i][i];        
        }         
         return newresult;    
          }
	
	public double[][] solve(double[][] a, double[][] b) {
		int nm = a.length;
		int n=a[0].length;
		int nn=b[0].length;
		double[][] x = new double[n][nn];
		for(int i=0;i<nn;i++){
			double[] bb = new double[n];
			for(int j=0;j<nm;j++) bb[j] = b[j][i];
			double[][] xishu = new double[nm][n+1]; 
			for(int p=0;p<nm;p++){
				xishu[p][n] = bb[p];
				for(int q=0;q<n;q++) xishu[p][q] = a[p][q];
			}		        
            double[] res = caculate(xishu);
            for(int p=0;p<n;p++) x[p][i] = res[p];
		}
            //double[][] xishu = { {4, 1, 2}, {2, 3, 5}};          
            //Test2 c = new Test2();         
            //double[] res = c.caculate(xishu);          
            //DecimalFormat df = new DecimalFormat("#.00");         
         return x;                      
        }
	
	public double[] backsolve(double[][] a, double[] b, boolean transpose) {
		int nm = a.length;
		int n=a[0].length;
		for(int i=0;i<nm;i++){
			for(int j=0;j<i;j++) a[i][j] = 0;
		}
		if(transpose){
			a=transpose(a);
		}
			double[][] xishu = new double[nm][n+1]; 
			for(int p=0;p<nm;p++){
				xishu[p][n] = b[p];
				for(int q=0;q<n;q++) xishu[p][q] = a[p][q];
			}		        
            double[] res = caculate(xishu);

            //double[][] xishu = { {4, 1, 2}, {2, 3, 5}};          
            //Test2 c = new Test2();         
            //double[] res = c.caculate(xishu);          
            //DecimalFormat df = new DecimalFormat("#.00");         
         return res;                      
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
	/*
	public static void main(String[] args) {
		double [][] a = {{2.4,0.8,1.4},{0,2.2,1.0},{0,0,2.3}};
		double[] b = {8,-4,2};
		double[] x = backsolve(a,b,true);
		System.out.println(x[0]);
		System.out.println(x[1]);
		System.out.println(x[2]);
	}
	*/
}
