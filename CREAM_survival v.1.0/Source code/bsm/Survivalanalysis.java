package bsm;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.swing.JTextArea;

import bsm.core.Util;
import javastat.survival.inference.LogRankTest;
import javastat.util.DataManager;

public class Survivalanalysis {
	List<Module> ModuleSet;
	
	public Survivalanalysis(BSM_DataSet expressiondata, RegulatorBindingData pathwayData, 
			HashMap<String, List<Integer>> resultlist, HashMap<String, List<Integer>> posilist, 
			HashMap<String, List<Integer>> negalist, int groupcount, double[][] survivalfile, 
			JTextArea casetext) throws Exception {
		
		casetext.append(" Gene: "+expressiondata.numrows+"\t"+" Patients: "+expressiondata.numcols+"\n");
		casetext.append(" Running..."+"\n");
	    casetext.paintImmediately(casetext.getBounds());
	    
		ModuleSet = new ArrayList<Module>();
		EnrichmentAnalysis EA = new EnrichmentAnalysis();
		 
		resultlist.forEach((moduleid, list) -> {
			Module module = new Module();
			module.id = moduleid;
			module.genelist = list;
			module.posilist = posilist.get(moduleid);
			module.negalist = negalist.get(moduleid);
			
			casetext.append(" Processing module: "+module.id+"\t"+module.posilist.size()+ " positively correlated genes"+
			"\t"+module.negalist.size()+" negatively correlated genes"+"\n");
			casetext.paintImmediately(casetext.getBounds());
			
			 if(pathwayData != null && pathwayData.reg2GeneDedupliIndex.length > 0){
				 module.pathwayenrichment = EA.enrichment(module.posilist, module.negalist,
						 pathwayData.regNames, pathwayData.reg2GeneDedupliIndex, expressiondata.numrows);
			 }
			
			double[][] exp = new double[list.size()][expressiondata.numcols];
			for(int i=0;i<module.posilist.size();i++) {
				exp[i] = expressiondata.controlnorm[module.posilist.get(i)];
			}
			for(int i=0;i<module.negalist.size();i++) {
				exp[module.posilist.size()+i] = expressiondata.controlnorm[module.negalist.get(i)];
			}
			String[][] avgstd = new String[expressiondata.numcols][2];
			for(int i=0;i<exp[0].length;i++) {
				avgstd[i][0] = i+"";
				double Sum = 0;
				double Sumsq = 0;
				for(int j=0;j<module.posilist.size();j++) {
					Sum += exp[j][i];
					Sumsq += Math.pow(exp[j][i], 2);
				}
				for(int j=module.posilist.size();j<list.size();j++) {
					Sum += -exp[j][i];
					Sumsq += Math.pow(exp[j][i], 2);
				}
				double mean = Sum/list.size();
				avgstd[i][1] = Math.sqrt(Sumsq/list.size() - Math.pow(mean,2))+"";
			}
	
			module.stdlist = Util.BubbleSort_inc(avgstd, avgstd.length, 1);
			
			for(int i=0;i<groupcount;i++) {
				module.highsamples.add(Integer.parseInt(module.stdlist[i][0]));
				module.lowsamples.add(Integer.parseInt(module.stdlist[avgstd.length-i-1][0]));
			}
			double[] timehigh = new double[groupcount];
			double[] timelow = new double[groupcount];
			double[] censorhigh = new double[groupcount];
			double[] censorlow = new double[groupcount];
			for(int i=0;i<groupcount;i++) {
				censorhigh[i] = survivalfile[module.highsamples.get(i)][0];
				timehigh[i] = survivalfile[module.highsamples.get(i)][1];
				censorlow[i] = survivalfile[module.lowsamples.get(i)][0];
				timelow[i] = survivalfile[module.lowsamples.get(i)][1];
			}
			DataManager dm = new DataManager();
	        LogRankTest testclass1 = new LogRankTest(timehigh, censorhigh, timelow, censorlow);
	        module.pvalue = dm.roundDigits(testclass1.pValue, 3.0);
	        //System.out.println(module.id+"\t"+module.pvalue);
			
			if(module.posilist.size() > 0) {
				module.positable = new String[module.posilist.size()][groupcount*2+1];
				for(int i=0;i<module.posilist.size();i++) {
					module.positable[i][0] = expressiondata.genenames[module.posilist.get(i)];
					for(int j=0;j<groupcount;j++) {
						module.positable[i][j+1] = expressiondata.controldata[module.posilist.get(i)]
								[module.highsamples.get(j)]+"";
					}
					for(int j=0;j<groupcount;j++) {
						module.positable[i][j+groupcount+1] = expressiondata.controldata[module.posilist.get(i)]
								[module.lowsamples.get(j)]+"";
					}
				}
			}
			if(module.negalist.size() > 0) {
				module.negatable = new String[module.negalist.size()][groupcount*2+1];
				for(int i=0;i<module.negalist.size();i++) {
					module.negatable[i][0] = expressiondata.genenames[module.negalist.get(i)];
					for(int j=0;j<groupcount;j++) {
						module.negatable[i][j+1] = expressiondata.controldata[module.negalist.get(i)]
								[module.highsamples.get(j)]+"";
					}
					for(int j=0;j<groupcount;j++) {
						module.negatable[i][j+groupcount+1] = expressiondata.controldata[module.negalist.get(i)]
								[module.lowsamples.get(j)]+"";
					}
				}
			}
			ModuleSet.add(module);
        });
		
	}
	
	public class Module {
		String id;
		double pvalue;
		List<Integer> genelist;
		List<Integer> posilist;
		List<Integer> negalist;
		List<Integer> highsamples;
		List<Integer> lowsamples;
		String[][] positable;
		String[][] negatable;
		String[][] stdlist;
		String[][] tfenrichment;
		String[][] pathwayenrichment;
		String[][] mirnaenrichment;
		
		Module() {
			initmodule();
		}
		void initmodule() {
			genelist = new ArrayList<Integer>();
			posilist = new ArrayList<Integer>();
			negalist = new ArrayList<Integer>();
			highsamples = new ArrayList<Integer>();
			lowsamples = new ArrayList<Integer>();
			//stdlist = new String[expressiondata.numcols][2];
		}
	}

}
