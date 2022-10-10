package ssclust;
/**
 * ����������ת��Ϊ�����Ǿ����һ�����������磺A=QR
 * */
public class ToDiagMatrix {
      private double[][] R;//�����Ǿ���
      private double[][] Q;//��������
      private double[][] A;//Ҫ�ֽ�ľ���A=QR
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
       * ��R�Խ����ϵ�Ԫ��
       * */
      public double Uii(int i){
    	  double sum=0.0;
    	  for(int j=0;j<n;j++){
    		  sum+=A[j][i]*A[j][i];
    	  }
    	  return Math.sqrt(sum);
      }
      /**��r���ǶԽ����ϵ�Ԫ��*/
      public double Uij(int i,int j){
    	  double u=0.0;
    	  for(int k=0;k<n;k++){
    		  u+=Q[k][i]*A[k][j];
    	  }
    	  return u;
      }
      /**�����Q��i������*/
      public void setQi(int i){
    	  for(int j=0;j<n;j++){
    		  Q[j][i]=A[j][i]/R[i][i];
    	  }
      }
      /**��I�θ���A*/
      public void updateA(int i){//��a1=a1-u01*q0
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
