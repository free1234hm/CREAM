package ssclust;

import java.util.Scanner;

//三样条插值
public class Spline {

	double h[],u[],q[],g[],M[];
	double S[];
	double x[],y[];	//插值节点、对应的节点值
	int n;
	//----------------追赶法法中需要的参数-----------
	double b[];//系数矩阵组成b,
	//系数
	double r[],a1[],b1[];
	//求解
	double x0[],y0[];
	public Spline(double []x,double []y,int n)
	{
		this.x=x;
		this.y=y;
		this.n=n;
	}
	//追赶法
	public double[] zhuigan(double a[],double c[],double d[])
	{  
		//Ax=d;d方程右边的矩阵,A为系数矩阵	
		//a系数矩阵组成a,c系数矩阵组成c,n节点数
		b=new double[n];//系数矩阵组成b,
	    //System.out.println(a.length+"   "+c.length);
		//系数
		r=new double[n];
		a1=new double[n];
		b1=new double[n];
		//求解
		 x0=new double[n];
		y0=new double[n];
		//赋值a/b/c
		for(int i=0;i<n;i++)
		{
			b[i]=2;
		}
		
		//三角分解
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
		 //求解方程Ax=b;
		 //1.Ly=b;
		 y0[0]=d[0]/a1[0];
		 for(int i=1;i<n-2;i++)
		 {
			 y0[i]=(d[i]-r[i]*y0[i-1])/a1[i];
			 //System.out.println("M的计算中出现问题了吗？?? "+y0[i]);
		 }
		 //2.Ux=y;
		 x0[n-2]=y0[n-2];
		 for (int i=n-3;i>=0;i--)
		 {
			 x0[i]=y0[i]-b1[i]*x0[i+1];
		 }
		 return x0;
	}
	//初始化h,u,q
	public double[] a( )
	{ //x，y,n分别为节点、节点对应的函数值、节点数
	   h=new double[n-1];
	   u=new double[n-2];
	   q=new double[n-2];
	   g=new double[n];
		//x为节点，y为对应节点x的节点值,n节点长度
		//AM=G
		for(int i=0;i<n-1;i++)
		{
			h[i]=x[i+1]-x[i];//x的导数
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
	 M=zhuigan(u,q,g);//调用追赶法计算M
	 return M;
	}
	//通过插值函数S计算对应yy中的节点的节点值
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

