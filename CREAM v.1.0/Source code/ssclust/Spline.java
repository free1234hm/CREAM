package ssclust;

import java.util.Scanner;

//��������ֵ
public class Spline {

	double h[],u[],q[],g[],M[];
	double S[];
	double x[],y[];	//��ֵ�ڵ㡢��Ӧ�Ľڵ�ֵ
	int n;
	//----------------׷�Ϸ�������Ҫ�Ĳ���-----------
	double b[];//ϵ���������b,
	//ϵ��
	double r[],a1[],b1[];
	//���
	double x0[],y0[];
	public Spline(double []x,double []y,int n)
	{
		this.x=x;
		this.y=y;
		this.n=n;
	}
	//׷�Ϸ�
	public double[] zhuigan(double a[],double c[],double d[])
	{  
		//Ax=d;d�����ұߵľ���,AΪϵ������	
		//aϵ���������a,cϵ���������c,n�ڵ���
		b=new double[n];//ϵ���������b,
	    //System.out.println(a.length+"   "+c.length);
		//ϵ��
		r=new double[n];
		a1=new double[n];
		b1=new double[n];
		//���
		 x0=new double[n];
		y0=new double[n];
		//��ֵa/b/c
		for(int i=0;i<n;i++)
		{
			b[i]=2;
		}
		
		//���Ƿֽ�
		for(int i=0;i<n-2;i++)
		{
			r[i]=a[i];
			
		}
		 a1[0]=b[0];
		 b1[0]=c[0]/a1[0];
		 for(int i=1;i<n-2;i++)
		 { 
			 a1[i]=b[i]-r[i]*b1[i-1];
			 b1[i]=c[i]/a1[i];
		 }
		 //��ⷽ��Ax=b;
		 //1.Ly=b;
		 y0[0]=d[0]/a1[0];
		 for(int i=1;i<n-2;i++)
		 {
			 y0[i]=(d[i]-r[i]*y0[i-1])/a1[i];
			 //System.out.println("M�ļ����г�����������?? "+y0[i]);
		 }
		 //2.Ux=y;
		 x0[n-2]=y0[n-2];
		 for (int i=n-3;i>=0;i--)
		 {
			 x0[i]=y0[i]-b1[i]*x0[i+1];
		 }
		 return x0;
	}
	//��ʼ��h,u,q
	public double[] a( )
	{ //x��y,n�ֱ�Ϊ�ڵ㡢�ڵ��Ӧ�ĺ���ֵ���ڵ���
	   h=new double[n-1];
	   u=new double[n-2];
	   q=new double[n-2];
	   g=new double[n];
		//xΪ�ڵ㣬yΪ��Ӧ�ڵ�x�Ľڵ�ֵ,n�ڵ㳤��
		//AM=G
		for(int i=0;i<n-1;i++)
		{
			h[i]=x[i+1]-x[i];//x�ĵ���
		}
		for(int i=0;i<n-2;i++)
		{
			u[i]=h[i]/(h[i] + h[i+1]);
		}
		for(int i=0;i<n-2;i++)
		{
			q[i]=1-u[i];
		}
		double y0=(y[1]-y[0])/(x[1]-x[0]);
		double y3=(y[n-1]-y[n-2])/(x[n-1]-x[n-2]);
		g[0]=6/h[0]*((y[1]-y[0])/h[0]-y0);
		for(int i=1;i<n-1;i++)
		{
			g[i]=6/(h[i-1]+h[i])*((y[i+1]-y[i])/h[i]-((y[i]-y[i-1])/h[i-1]));
		}
		g[n-1]=6/h[n-2]*(y3-(y[n-1]-y[n-2])/h[n-2]);
	 M=zhuigan(u,q,g);//����׷�Ϸ�����M
	 return M;
	}
	//ͨ����ֵ����S�����Ӧyy�еĽڵ�Ľڵ�ֵ
	public double[] s(int []yy)
	{  
		double M[]=a();
		S=new double[yy[yy.length-1]+1];
	
			for(int i=0;i<n-1;i++)
			{ 
				for(int f=0;f<yy.length;f++)
			{
				if(yy[f]>x[i]&&yy[f]<x[i+1])
				{ 
					S[yy[f]]=Math.pow((x[i+1]-yy[f]),3)*M[i]/(6*h[i])+Math.pow((yy[f]-x[i]),3)*M[i+1]/(6*h[i]) + (x[i+1]-yy[f])*(y[i]-h[i]*h[i]*M[i]/6)/h[i] + (yy[f]-x[i])*(y[i+1]-h[i]*h[i]*M[i+1]/6)/h[i];
					//System.out.println("S[" +f+"]="+S);
				}
				else if(yy[f]==x[i])
				{
					S[yy[f]]=y[i];
				//	i++;
				}
				else if(yy[f]==x[i+1])
				{
					S[yy[f]]=y[i+1];
					//i++;
				}
			}
			}return S;		
	}	
	}

