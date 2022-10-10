package bsm;

import java.util.Random;

import ssclust.Pearsonr;

public class Corrtest {
	public static void main(String[] args) throws Exception {
		Pearsonr pr = new Pearsonr();
		Random r = new Random();
    	double[][] ran = {{1,2,3},
    			          {3,2,2},
    			          {4,6,3},
    			          {5,7,6},
    			          {6,6,6},
    			          {1,2,3},
    			          {3,2,2},
    			          {4,6,3},
    			          {5,7,6},
    			          {6,6,6},
    			          {1,2,3},
    			          {3,2,2},
    			          {7,6,5},
    			          {2,7,3}};
    	
    	double[][] rr = pr.pearsonp2(ran);
    	double[][] pp = pr.pearsonp3(ran);
    	
    	for(int i=0;i<rr.length;i++){
    		for(int j=0;j<rr[0].length;j++){
    			System.out.print(rr[i][j]+"\t");
    		}
    		System.out.println();
    	}
    	System.out.println("===================================================================");
    	for(int i=0;i<pp.length;i++){
    		for(int j=0;j<pp[0].length;j++){
    			System.out.print(pp[i][j]+"\t");
    		}
    		System.out.println();
    	}
    

		
	}
}
