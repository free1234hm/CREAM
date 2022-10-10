package ssclust;


public class Matrix_dot {

public double[][] dot(double[][] matrix1, double[][] matrix2) {

if(matrix1.length!=matrix2[0].length){//若无法相乘则退出
System.out.println("ivalid input");
System.exit(0);
}

double[][] r = new double[matrix1[0].length][matrix2.length];
for(int i=0;i<r.length;++i){
for(int j=0;j<r[i].length;++j){//每一个r[i][j]的运算：
r[i][j]=0;//初始化
for(int k=0;k<matrix2.length;++k)
r[i][j]+=matrix1[i][k]*matrix2[k][j];
}
}
//输出结果
return r;
}
}
