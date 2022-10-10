package bsm;


import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.swing.JTextArea;

import bsm.core.DC_pvalue;

public class DCanalysis {
	List<Module> ModuleSet;
	
	public DCanalysis(BSM_DataSet controldata, BSM_DataSet casedata, RegulatorBindingData pathwayData, 
			HashMap<String, List<Integer>> resultlist, HashMap<String, List<Integer>> posilist, 
			HashMap<String, List<Integer>> negalist, JTextArea casetext) throws Exception {
		
		casetext.append(" Gene: "+controldata.numrows+"\t"+" Control samples: "+controldata.numcols
				+" Case samples: "+casedata.numcols+"\n");
		casetext.append(" Running..."+"\n");
	    casetext.paintImmediately(casetext.getBounds());
	    
		ModuleSet = new ArrayList<Module>();
		EnrichmentAnalysis EA = new EnrichmentAnalysis();
		
		resultlist.forEach((moduleid, list) -> {
			Module module = new Module();
			module.id = moduleid;
			module.list = list;
			module.posilist = posilist.get(moduleid);
			module.negalist = negalist.get(moduleid);
			
			if(module.posilist == null) {
				module.posilist = new ArrayList<Integer>();
			}
			if(module.negalist == null) {
				module.negalist = new ArrayList<Integer>();
			}
			
			casetext.append(" Analyzing module: "+module.id+"\t"+module.posilist.size()+ " positively correlated genes and "+
			+module.negalist.size()+" negatively correlated genes"+"\n");
			casetext.paintImmediately(casetext.getBounds());
			
			 if(pathwayData != null && pathwayData.reg2GeneDedupliIndex.length > 0){
				 module.pathwayenrichment = EA.enrichment(module.posilist, module.negalist,
						 pathwayData.regNames, pathwayData.reg2GeneDedupliIndex, controldata.numrows);
			 }
			 int samplenum = Math.min(controldata.numcols, casedata.numcols);
			 
			 DC_pvalue dc = new DC_pvalue();
			 module.dcpvalue = dc.gcpathway(controldata, casedata, list, samplenum);
			 //module.dcpvalue = new double[3];
			 module.controltable = new String[module.list.size()][controldata.numcols+1];
			 module.casetable = new String[module.list.size()][casedata.numcols+1];
			 for(int i=0;i<module.posilist.size();i++) {
				 module.controltable[i][0] = controldata.genenames[module.posilist.get(i)];
				 module.casetable[i][0] = casedata.genenames[module.posilist.get(i)];
				 for(int j=0;j<controldata.numcols;j++) {
					 module.controltable[i][j+1] = controldata.controldata[module.posilist.get(i)][j]+"";
				 }
				 for(int j=0;j<casedata.numcols;j++) {
					 module.casetable[i][j+1] = casedata.controldata[module.posilist.get(i)][j]+"";
				 }
			 }
			 for(int i=0;i<module.negalist.size();i++) {
				 module.controltable[i+module.posilist.size()][0] = controldata.genenames[module.negalist.get(i)];
				 module.casetable[i+module.posilist.size()][0] = casedata.genenames[module.negalist.get(i)];
				 for(int j=0;j<controldata.numcols;j++) {
					 module.controltable[i+module.posilist.size()][j+1] = controldata.controldata[module.negalist.get(i)][j]+"";
				 }
				 for(int j=0;j<casedata.numcols;j++) {
					 module.casetable[i+module.posilist.size()][j+1] = casedata.controldata[module.negalist.get(i)][j]+"";
				 }
			 }
			ModuleSet.add(module);
        });
	}
	
	public class Module {
		String id;
		double pvalue;
		List<Integer> list;
		List<Integer> posilist;
		List<Integer> negalist;
		String[][] controltable;
		String[][] casetable;
		double[] dcpvalue;
		String[][] pathwayenrichment;
		
		Module() {
			initmodule();
		}
		void initmodule() {
			list = new ArrayList<Integer>();
			posilist = new ArrayList<Integer>();
			negalist = new ArrayList<Integer>();
		}
	}

}
