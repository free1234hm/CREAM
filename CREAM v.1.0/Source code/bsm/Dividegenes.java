package bsm;



import java.util.ArrayList;
import java.util.List;
import ssclust.Pearsonr;
import bsm.core.Util;



public class Dividegenes {
	
	public List<Integer> dividegene(String[][] express, double[] meancurve, 
			double minicor, Integer minnum){
		
    	int row = express.length;
    	int col = express[0].length-1;
    	List<Integer> list = new ArrayList<Integer>();
    	List<Double> pearsonlist = new ArrayList<Double>();
		
		for(int i=0;i<row;i++){
			int geneid = Integer.parseInt(express[i][0]);
			Double[] value = new Double[col];
			for(int j=0;j<col;j++) value[j] = Double.parseDouble(express[i][j+1]);
			
			Pearsonr pearson = new Pearsonr();
			double r = pearson.pearsonr(value, meancurve);	
			
			if(Math.abs(r)<minicor) {
				list.add(geneid);
				pearsonlist.add(Math.abs(r));
			}
		}
		
	if(list.size() >= minnum){
		 if((row-list.size()) >= minnum){
			 return list;
		 }else{
			 String[][] value = new String[list.size()][2];
			 for(int i=0;i<value.length;i++){
				 value[i][0] = list.get(i)+"";
				 value[i][1] = pearsonlist.get(i)+"";
			 }
			 Util.BubbleSort_inc(value, value.length, 1);
			 List<Integer> list2 = new ArrayList<Integer>();
			 int len = value.length/2;
			 for(int i=0;i<len;i++){
				 list2.add(Integer.parseInt(value[i][0]));
			 }
			 return list2;
		 }	
	  }else{
		return null;
	  }
	}
	
}

