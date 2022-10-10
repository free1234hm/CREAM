package ssclust;
/**
 * 把正定矩阵转化为上三角矩阵和一个正交矩阵，如：A=QR
 * */
public class ToDiagMatrix {
      private double[][] R;//上三角矩阵
      private double[][] Q;//正交矩阵
      private double[][] A;//要分解的矩阵A=QR
      private int n;
      public ToDiagMatrix(double[][] A){
    	  this.A=A;
    	  n=A.length;
          this.R=new double[n][n];
          this.Q=new double[n][n];
          this.initR();
          schmidt();
      }
      public void initR(){
    	  for(int i=0;i<n;i++){
    		  for(int j=0;j<n;j++){
    			  R[i][j]=0.0;
    		  }
    	  }
      }
      public void schmidt(){
    	  for(int i=0;i<n;i++){
    		  R[i][i]=this.Uii(i);
    		  this.setQi(i);
    		  for(int j=i+1;j<n;j++){
    			  R[i][j]=this.Uij(i, j);
    		  }
    		  this.updateA(i);
    	  }
      }
      /**
       * 求R对角线上的元素
       * */
      public double Uii(int i){
    	  double sum=0.0;
    	  for(int j=0;j<n;j++){
    		  sum+=A[j][i]*A[j][i];
    	  }
    	  return Math.sqrt(sum);
      }
      /**求r不是对角线上的元素*/
      public double Uij(int i,int j){
    	  double u=0.0;
    	  for(int k=0;k<n;k++){
    		  u+=Q[k][i]*A[k][j];
    	  }
    	  return u;
      }
      /**求矩阵Q的i列向量*/
      public void setQi(int i){
    	  for(int j=0;j<n;j++){
    		  Q[j][i]=A[j][i]/R[i][i];
    	  }
      }
      /**第I次更新A*/
      public void updateA(int i){//如a1=a1-u01*q0
    	  for(int j=i+1;j<n;j++){
    		  for(int k=0;k<n;k++){
    			  A[k][j]=A[k][j]-R[i][j]*Q[k][i];
    		  }
    	  }
      }
	public double[][] getR() {
		return R;
	}
	public double[][] getQ() {
		return Q;
	}
}
