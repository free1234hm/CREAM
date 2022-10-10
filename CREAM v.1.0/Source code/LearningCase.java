package bsm;

import bsm.core.*;
import javax.swing.*;
import org.apache.commons.math3.distribution.NormalDistribution;
import ssclust.BuildSSM;
import ssclust.Pearsonr;
import java.util.*;
import java.util.List;
import java.text.*;
import java.io.*;

/**
 * This class implements the core methods for learning the BSM models
 */

public class LearningCase {

	static final boolean BDEBUG = false;
	static final boolean BDEBUGMODEL = false;
	boolean bhasmerge = false;
	int MinNum = 10;
	int minnum = 5;
	int npath;
	int nrandomseed;
	int ninitsearchval;
	int MinGeneNum;
	JButton currentButton;
	JTextArea casetext;
	static final String SZDELIM = "|;,";
	double[] dprevouterbestlog;
	double BEPSILON = 0.00001;
	double EPSILON = 0.005;
	int nmaxchild;
	int nminchild;
	double trainPearson;
	double minPearson = 0.4;
	int numcols;
	int numrows;
	int MINPATH = 5;
	
	static final double RIDGE = 1;
	/**
	 * Object containing all of the Regulator-Gene binding data.
	 */
	List<RegulatorBindingData> pathwayData;
	double[][] traindata;
	boolean busego;
	double NDmean;
	double NDsigma;
	NormalDistribution Normalcorr;
	double Sumcorr;   
	double Sumcorrsq;
	Casenode treeptr;
	NumberFormat nf4;
	double dbestlog;
	double peaklike;
	double currentbestlike;
	Casenode bestTree = null;
	BSM_DataSet theDataSet; 
	int numtotalPath = 1;
	List<List<String>> savedmodel;
	String[][] resultList;
	List<String[][]> minenrichmenttable;
	boolean[] IsSeed;

	public LearningCase(BSM_DataSet theDS, List<RegulatorBindingData> pathwayData,
			boolean busego, List<List<String>> savedmodel, String szmaxchild, String szminchild,
			String Pearson1, String epsilon, JTextArea casetext, JButton endSearchButton) throws Exception {
	
		nf4 = NumberFormat.getInstance(Locale.ENGLISH); //nf3 transform double to 4 decimal places
		nf4.setMinimumFractionDigits(6);
		nf4.setMaximumFractionDigits(6);
		this.theDataSet = theDS;   //The time series data
		this.nmaxchild = Integer.parseInt(szmaxchild);
		this.nminchild = Integer.parseInt(szminchild);
		this.pathwayData = pathwayData;
		this.busego = busego;
		this.casetext = casetext;
		this.savedmodel = savedmodel;
		this.numrows = theDataSet.numrows;
		this.numcols = theDataSet.numcols;
		this.trainPearson = Double.parseDouble(Pearson1);
		this.EPSILON = Double.parseDouble(epsilon);
		treeptr = new Casenode(); //Node Initialization
		
	    casetext.append(" Gene: "+numrows+"\t"+" Sample: "+numcols+"\n");
	    casetext.paintImmediately(casetext.getBounds());
		
		traindata = theDataSet.controlnorm;
		
	    if(savedmodel != null && savedmodel.size()>0){
	    	System.out.println("Module number: "+savedmodel.size());
		    generateTree(treeptr, savedmodel);
		    trainhmm2(treeptr);
		    EnrichmentAnalysis EA = new EnrichmentAnalysis();

		    softclustering(treeptr);
		    
		    if(pathwayData != null && pathwayData.size() > 0){
		    	for(int i=0;i<pathwayData.size();i++) {
		    		RegulatorBindingData knownmodule = pathwayData.get(i);
		    		for(int nchild=0;nchild<treeptr.numchildren;nchild++){
			    		List<Integer> posilist = treeptr.nextptr[nchild].posi_genelist;
						List<Integer> negalist = treeptr.nextptr[nchild].nega_genelist;
						treeptr.nextptr[nchild].enrichment.add(EA.enrichment(posilist, negalist,
								knownmodule.regNames, knownmodule.reg2GeneDedupliIndex, numrows));
			    	}
		    	}
		    }
		    
		    try{
				BufferedWriter outXml = new BufferedWriter(new FileWriter("Final_modules.txt"));
				outXml.write("Module"+"\t"+"Gene"+"\t"+"Correlation"+"\t");
				for(int i=0;i<theDataSet.numcols-1;i++) {
					outXml.write(theDataSet.dsamplemins[i]+"\t");
				}
				outXml.write(theDataSet.dsamplemins[theDataSet.numcols-1]+"\n");
				
				for(int i=0;i<treeptr.numchildren;i++){
					List<Integer> posigenelist = treeptr.nextptr[i].posi_genelist;
					List<Integer> negagenelist = treeptr.nextptr[i].nega_genelist;
					if(posigenelist.size()>0){
						for(int j=0;j<posigenelist.size();j++){
							int geneindex = posigenelist.get(j);
					    	outXml.write("Module_"+(i+1)+"\t"+theDataSet.genenames[geneindex]+"\t"+"1"+"\t");
					    	for(int m=0;m<theDataSet.numcols-1;m++) {
					    		outXml.write(theDataSet.controldata[geneindex][m]+"\t");
					    	}
					    	outXml.write(theDataSet.controldata[geneindex][theDataSet.numcols-1]+"\n");
					    }
					}
					if(negagenelist.size()>0){
						for(int j=0;j<negagenelist.size();j++){
							int geneindex = negagenelist.get(j);
					    	outXml.write("Module_"+(i+1)+"\t"+theDataSet.genenames[geneindex]+"\t"+"-1"+"\t");
					    	for(int m=0;m<theDataSet.numcols-1;m++) {
					    		outXml.write(theDataSet.controldata[geneindex][m]+"\t");
					    	}
					    	outXml.write(theDataSet.controldata[geneindex][theDataSet.numcols-1]+"\n");
					    }
					}
				}
				outXml.flush(); 
				outXml.close();
				System.out.println("DONE");
			}catch (Exception e) {
				System.out.println("FALSE"); 
			e.printStackTrace(); 
			}

		    
	    }else{
	    	List<Integer> genelist = new ArrayList<Integer>();
		    for(int i=0;i<numrows;i++) genelist.add(i);
		    BuildSSM ssm = new BuildSSM();
		    double[] curve = ssm.meancurve(traindata, genelist);
		    buildEmptyTree(treeptr, curve);
		    endSearchButton.setEnabled(true);
		    
		    searchstage1();
		    trainhmm2(treeptr);

		    if (endSearchButton != null) {
				endSearchButton.setEnabled(false);
			}
			casetext.append("================================"+"\n");
			casetext.paintImmediately(casetext.getBounds());
			
			System.out.println("Final path = " + treeptr.numchildren);
			try{
				BufferedWriter outXml = new BufferedWriter(new FileWriter("D:/finalmodel.txt"));
				for(int i=0;i<treeptr.numchildren;i++){
					 double[] mean = treeptr.nextptr[i].curve;
					for(int j=0;j<mean.length-1;j++){
				    	outXml.write(mean[j] + "\t");
				    }
					outXml.write(mean[mean.length-1]+"\n");
				}
				outXml.flush(); 
				outXml.close();
				System.out.println("DONE");
				
			} catch (Exception e) {
				System.out.println("FALSE"); 
			    e.printStackTrace(); 
			}
			
			EnrichmentAnalysis EA = new EnrichmentAnalysis();
		    
		    softclustering(treeptr);
		    
		    if(pathwayData != null && pathwayData.size() > 0){
		    	for(int i=0;i<pathwayData.size();i++) {
		    		RegulatorBindingData knownmodule = pathwayData.get(i);
		    		for(int nchild=0;nchild<treeptr.numchildren;nchild++){
			    		List<Integer> posilist = treeptr.nextptr[nchild].posi_genelist;
						List<Integer> negalist = treeptr.nextptr[nchild].nega_genelist;
						treeptr.nextptr[nchild].enrichment.add(EA.enrichment(posilist, negalist,
								knownmodule.regNames, knownmodule.reg2GeneDedupliIndex, numrows));
			    	}
		    	}
		    }
		    
		    try{
				BufferedWriter outXml = new BufferedWriter(new FileWriter("Final_modules.txt"));
				outXml.write("Module"+"\t"+"Gene"+"\t"+"Correlation"+"\t");
				for(int i=0;i<theDataSet.numcols-1;i++) {
					outXml.write(theDataSet.dsamplemins[i]+"\t");
				}
				outXml.write(theDataSet.dsamplemins[theDataSet.numcols-1]+"\n");
				
				for(int i=0;i<treeptr.numchildren;i++){
					List<Integer> posigenelist = treeptr.nextptr[i].posi_genelist;
					List<Integer> negagenelist = treeptr.nextptr[i].nega_genelist;
					if(posigenelist.size()>0){
						for(int j=0;j<posigenelist.size();j++){
							int geneindex = posigenelist.get(j);
					    	outXml.write("Module_"+(i+1)+"\t"+theDataSet.genenames[geneindex]+"\t"+"1"+"\t");
					    	for(int m=0;m<theDataSet.numcols-1;m++) {
					    		outXml.write(theDataSet.controldata[geneindex][m]+"\t");
					    	}
					    	outXml.write(theDataSet.controldata[geneindex][theDataSet.numcols-1]+"\n");
					    }
					}
					if(negagenelist.size()>0){
						for(int j=0;j<negagenelist.size();j++){
							int geneindex = negagenelist.get(j);
					    	outXml.write("Module_"+(i+1)+"\t"+theDataSet.genenames[geneindex]+"\t"+"-1"+"\t");
					    	for(int m=0;m<theDataSet.numcols-1;m++) {
					    		outXml.write(theDataSet.controldata[geneindex][m]+"\t");
					    	}
					    	outXml.write(theDataSet.controldata[geneindex][theDataSet.numcols-1]+"\n");
					    }
					}
				}
				outXml.flush(); 
				outXml.close();
				System.out.println("DONE");
			}catch (Exception e) {
				System.out.println("FALSE"); 
			e.printStackTrace(); 
			}
	    }
	    
	    resultList = new String[treeptr.numchildren][numcols+4];
	    for(int i=0;i<treeptr.numchildren;i++){
			resultList[i][0] = i+1+"";
			resultList[i][1] = (treeptr.nextptr[i].posi_genelist.size()+
					treeptr.nextptr[i].nega_genelist.size())+"";
			resultList[i][2] = treeptr.nextptr[i].posi_genelist.size()+"";
			resultList[i][3] = treeptr.nextptr[i].nega_genelist.size()+"";
			double[] mean = treeptr.nextptr[i].curve;
			for(int j=0;j<mean.length;j++){
				resultList[i][j+4] = mean[j]+"";
		    }
		}
	    
		boolean bagain; //Delete modules with lower than 5 genes
			do {
				bagain = false;
				softclustering(treeptr);
				int theMinPathRec = findMinPath(treeptr);		
				if (MinGeneNum < MINPATH) {//删除5个基因以下的路径
					casetext.append(" Module "+theMinPathRec+" are filtered for including lower than 5 genes."+"\n");
					casetext.paintImmediately(casetext.getBounds());
					deleteMinPath(theMinPathRec, treeptr); //删除包含最少基因的路径
					bagain = true;
					trainhmm(treeptr);
					if (BDEBUG) {
						System.out.println("after retrain");
					}
					numtotalPath--;
				}
			} while (bagain);
			casetext.append("================================"+"\n");
			casetext.paintImmediately(casetext.getBounds());	
	}

	private void generateTree(Casenode ptr, List<List<String>> savedmodel){
		ptr.numchildren = savedmodel.size();
		for(int i=0;i<savedmodel.size();i++){
			List<String> list = savedmodel.get(i);
			double[] curve = new double[list.size()];
			if(curve.length != theDataSet.numcols){
				throw new IllegalArgumentException("The column number of saved model should be "
						+ "equal to the case samples.");
			}else{
				for(int j=0;j<curve.length;j++) curve[j] = Double.parseDouble(list.get(j));
				ptr.nextptr[i] = new Casenode(ptr);
				ptr.nextptr[i].curve = curve;	
			}
		}
	}
	
	private int findMinPath(Casenode ptr) throws Exception {
	    int nmin = 0;
	    MinGeneNum = Integer.MAX_VALUE;
		for (int nindex = 0; nindex < ptr.numchildren; nindex++) {
			int genecount = ptr.nextptr[nindex].posi_genelist.size() 
					+ ptr.nextptr[nindex].nega_genelist.size();
			if (genecount < MinGeneNum) {
				nmin = nindex;
				MinGeneNum = genecount;
			}
		}
		return nmin;
     }

	/**
	 * First step of tree searching，add and delete path
	 */
	public void searchstage1() throws Exception {
		dbestlog = trainhmm(treeptr);
		dprevouterbestlog = new double[nmaxchild];
		dprevouterbestlog[numtotalPath-1] = dbestlog;
		boolean bendsearchlocal;  //定义变量记录是否结束搜索
		Casenode finalTree = (Casenode) treeptr.clone();
		peaklike = 0;
		boolean keepon;
		do {
			bestTree = treeptr;
			keepon = false;
			/************************add path****************************/
			double avgcor = hardclustering(treeptr);
			System.out.println("Average corrlation = " + avgcor);
			double traincor = Math.max(0.4, Math.min(avgcor, 0.8));
			
			currentbestlike = 0;
			
			for(int nchild = 0; nchild < treeptr.numchildren; nchild++){
				List<Integer> posigenes = treeptr.nextptr[nchild].posi_genelist;
				List<Integer> negagenes = treeptr.nextptr[nchild].nega_genelist;

				if(posigenes.size() >= MinNum){
					List<Integer> posilist = new ArrayList<Integer>();
					List<Integer> restlist = new ArrayList<Integer>();
					for(int gene:posigenes){
						double r = treeptr.nextptr[nchild].dcorrec[gene];
						if(Math.abs(r) < traincor){
							posilist.add(gene);
						}else{
							restlist.add(gene);
						}
					}
					//for(int gene:negagenes){
					//	restlist.add(gene);
					//}
					if(posilist.size() >= minnum && restlist.size() >= minnum){
						traverseandadd(treeptr, nchild, posilist, restlist);
					}
				}
				
				if(negagenes.size() >= MinNum){
					List<Integer> negalist = new ArrayList<Integer>();
					List<Integer> restlist = new ArrayList<Integer>();
					for(int gene:negagenes){
						double r = treeptr.nextptr[nchild].dcorrec[gene];
						if(Math.abs(r) < traincor){
							negalist.add(gene);
						}else{
							restlist.add(gene);
						}
					}
					//for(int gene:posigenes){
					//	restlist.add(gene);
					//}
					if(negalist.size() >= minnum && restlist.size() >= minnum){
						traverseandadd(treeptr, nchild, negalist, restlist);
					}
				}
			}
			
			/****************************************************/
			if(currentbestlike > 0){
				treeptr = bestTree;
				numtotalPath++;
				dprevouterbestlog[numtotalPath-1] = currentbestlike;
				List<Integer> similarchild = findsimilar(treeptr);
				if(similarchild.size()>0){
					boolean delete = traverseanddelete(similarchild, treeptr);
					if (delete){
						treeptr = bestTree;
						numtotalPath--;
						dprevouterbestlog[numtotalPath-1] = currentbestlike;
					}
				}
			} else {
				break;
			}
			/****************************************************/
			
			/****************************************************/
			if(numtotalPath > 1){
				if(dprevouterbestlog[numtotalPath-1] < dprevouterbestlog[numtotalPath-2]) {
					peaklike = dprevouterbestlog[numtotalPath-2];
				}else if((dprevouterbestlog[numtotalPath-1] - dprevouterbestlog[numtotalPath-2]) / 
						dprevouterbestlog[numtotalPath-1] > EPSILON && dprevouterbestlog[numtotalPath-1] > peaklike){
					finalTree = (Casenode) bestTree.clone();
				}
			}
			/****************************************************/
			
			bendsearchlocal = Main_interface.bendsearch; //bend searching
			System.out.println(" Current path: "+numtotalPath+"; Best liklihood: "
			        +dbestlog+ "; Current liklihood: " + currentbestlike + "; pre liklihood: "+ dprevouterbestlog[Math.max(numtotalPath-2, 0)]);
			casetext.append(" Current path: "+numtotalPath+"; Best liklihood: "+dbestlog+ "; Current liklihood: " + currentbestlike+"\n");
			casetext.paintImmediately(casetext.getBounds());
			
			for(int path=numtotalPath-1; path>Math.max(numtotalPath-4, 0); path--){
				if((dprevouterbestlog[path] - dprevouterbestlog[path-1]) / dprevouterbestlog[path] > EPSILON){
					keepon = true;
					break;
				}
			}
		} while(!bendsearchlocal && (((keepon || numtotalPath==1) && numtotalPath < nmaxchild) || (numtotalPath < nminchild)));
		treeptr = finalTree;
	}
	
	public List<Integer> findsimilar(Casenode treeptr) {
		List<Integer> similarchild = new ArrayList<Integer>();
		Set<Integer> set = new HashSet<Integer>();
		double[][] countsimilar = new double[treeptr.numchildren][2];
		for(int i=0;i<countsimilar.length;i++) countsimilar[i][0] = i;
		Pearsonr pr = new Pearsonr();
		for(int i=0;i<treeptr.numchildren;i++){
			for(int j=i+1;j<treeptr.numchildren;j++){
				double[] mean1 = treeptr.nextptr[i].curve;
				double[] mean2 = treeptr.nextptr[j].curve;
				double r = pr.cosineSimilarity(mean1, mean2);
				if(Math.abs(r) >= minPearson){
					countsimilar[i][1] += Math.abs(r);
					countsimilar[j][1] += Math.abs(r);
					set.add(i);
					set.add(j);
				}
			}
		}
		if(set.size()>10){
			countsimilar = Util.BubbleSort_dec(countsimilar, countsimilar.length, 1);
			for(int i=0;i<10;i++) similarchild.add((int)countsimilar[i][0]);
		}else{
			similarchild = new ArrayList<>(set);
		}
		return similarchild;
	}


	/**
	 * Builds the initial tree which is just a single chain with mean and
	 * standard deviation for each node being the global mean and standard
	 * deviation at the time point
	 */
	public void buildEmptyTree(Casenode currnode, double[] curve) {
			currnode.numchildren = 1;
			currnode.nextptr[0] = new Casenode(currnode);
			currnode.nextptr[0].curve = curve;
			currnode.nextptr[0].parent = currnode;
	}

	/**
	 * Deletes a child from the specified path on a cloned version of the root
	 * node
	**/
	public Casenode deletepath(int path, Casenode root) {
		Casenode treeroot = (Casenode) root.clone();
		Casenode ptr = treeroot;

		for (int nj = path + 1; nj < ptr.numchildren; nj++) {
			ptr.nextptr[nj - 1] = ptr.nextptr[nj];
		}
		ptr.numchildren--;
		return treeroot;
	}
	
	/**
	 * 添加一个路径，并验证likehood是否有改变
	 */
	private void traverseandadd(Casenode origroot, Integer nchild, List<Integer> list1, List<Integer> list2) 
			throws Exception {
		if (!Main_interface.bendsearch && (origroot != null) && (origroot.numchildren < nmaxchild)) {
			//System.out.println(numtotalPath+" "+origroot.numchildren);
				Casenode splittree = splitnode(origroot, nchild, list1, list2); 
				if (splittree != null) {  //If succeed to add path
					double dlog = trainhmm(splittree);	
					if (dlog > currentbestlike) {
						bestTree = splittree;
						currentbestlike = dlog;
					}
					if (dlog > dbestlog) {
						dbestlog = dlog;
					}
			    }
		}
	}

	/**
	 * 在一个节点后添加一个分岔节点
	 */
	private Casenode splitnode(Casenode origroot, Integer nchild, List<Integer> list1, List<Integer> list2) {
		
		Casenode treeroot = (Casenode) origroot.clone();
		Casenode ptr = treeroot;
		
		BuildSSM ssm1 = new BuildSSM();
		double[] curve1 = ssm1.meancurve(traindata, list1); //new path
		BuildSSM ssm2 = new BuildSSM();
		double[] curve2 = ssm2.meancurve(traindata, list2, ptr.nextptr[nchild].dcorrec); //old path
		
		ptr.nextptr[nchild].curve = curve2;
	    ptr.numchildren++;
		ptr.nextptr[ptr.numchildren - 1] = new Casenode(ptr);
		ptr.nextptr[ptr.numchildren - 1].curve = curve1;
		
		return treeroot;
	}


	/**
	 * Deletes the specificed path from the model starting from ndesiredlevel
	 */
	public void deleteMinPath(int path, Casenode root) {
		Casenode ptr = root;
		for (int nj = path + 1; nj < ptr.numchildren; nj++) {
			ptr.nextptr[nj - 1] = ptr.nextptr[nj];
		}
		ptr.numchildren--;
	}

	/**
	 * Helper function that searches for the best path to delete
	 */
	private boolean traverseanddelete(List<Integer> similarlist, Casenode origroot) throws Exception {
		double bestdeletlog = dprevouterbestlog[numtotalPath-2];
		boolean delete = false;
			for (int nindex:similarlist) {
				Casenode deletetree = deletepath(nindex, origroot);
				if (deletetree != null) {
					double dlog = trainhmm(deletetree);			
					if (dlog > bestdeletlog)
					{
						bestdeletlog = dlog;
						bestTree = deletetree;
						dbestlog = dlog;
						currentbestlike = dlog;
						delete = true;
					}
				}
			}
			return delete;
	}
	
	/**
	 * A Casenode corresponds to a state in the model
	 */
	public class Casenode {
		int ncurrtime;// use to count visited
		double[] curve;
		
		List<Integer> posi_genelist;
		List<Integer> nega_genelist; 
		double[] dcorrec;
		int numchildren;
		
		List<String[][]> enrichment;
		
		double[] dEsum;
		double[] dPsum;
		
		/** pointer to each of the next states */
		Casenode[] nextptr;
		/** pointer back to the parent */
		Casenode parent;
		
		/**
		 * Calls inittreenode with a null parent
		 */
		Casenode() {
			inittreenode(null);
		}

		Casenode(Casenode parent) {
			inittreenode(parent);
		}

		/**
		 * Calls inittreenode with parent
		 */
		void inittreenode(Casenode parent) {
			curve = new double[numcols];
			posi_genelist = new ArrayList<Integer>();
			nega_genelist = new ArrayList<Integer>();
			enrichment = new ArrayList<String[][]>();
			dEsum = new double[numcols];
			dPsum = new double[numcols];
			dcorrec = new double[numrows];
			this.parent = parent;
			nextptr = new Casenode[nmaxchild];
			numchildren = 0;
			for (int nindex = 0; nindex < nextptr.length; nindex++) {
				nextptr[nindex] = null;
			}
		}

		/**
		 * For making copy of nodes
		 */
		public Object clone() {
		    Casenode tnode = new Casenode();
			tnode.curve = new double[curve.length];
			for(int j=0;j<tnode.curve.length;j++) tnode.curve[j] = curve[j];
			tnode.dcorrec = new double[dcorrec.length];
			for(int j=0;j<tnode.dcorrec.length;j++) tnode.dcorrec[j] = dcorrec[j];
			tnode.numchildren = numchildren;
			tnode.parent = null;
			tnode.nextptr = new Casenode[nmaxchild]; //定义当前节点的子节点组
			for (int nindex = 0; nindex < nextptr.length; nindex++) {
				if (nextptr[nindex] == null) {
					tnode.nextptr[nindex] = null;
				} else {
					// cloning a new node
					tnode.nextptr[nindex] = (Casenode) nextptr[nindex].clone(); //null
					tnode.nextptr[nindex].parent = tnode;
				}
			}
			return tnode;  //返回当前节点tnode
		}
	}
	
	/**
	 * Clears any prior assignments of genes to paths through the model
	 */
	public void clearCounts(Casenode ptr){
		if (ptr != null) {
			ptr.posi_genelist.clear();
			ptr.nega_genelist.clear();
			for (int nchild = 0; nchild < ptr.numchildren; nchild++) {
				clearCounts(ptr.nextptr[nchild]);
			}
		}
	}
	
	public double hardclustering(Casenode treeptr){
		clearCounts(treeptr);
		double avgcorr = 0;
		for (int nrow = 0; nrow < numrows; nrow++) {
			double bestcorr = 0;
			int bestchild = -1;
			for (int nchild = 0; nchild < treeptr.numchildren; nchild++){
				double corr = treeptr.nextptr[nchild].dcorrec[nrow];
				if(Math.abs(corr) > Math.abs(bestcorr)){
					bestcorr = corr;
					bestchild = nchild;
				}
			}
			avgcorr += Math.abs(bestcorr);
			if(bestcorr > 0){
				treeptr.nextptr[bestchild].posi_genelist.add(nrow);
			}else if(bestcorr < 0){
				treeptr.nextptr[bestchild].nega_genelist.add(nrow);
			}
		}
		avgcorr = avgcorr/numrows;
		return avgcorr;
	}
	
	public void softclustering(Casenode treeptr){
		clearCounts(treeptr);
		IsSeed = new boolean[numrows];
		double[][] recweight = new double[numrows][treeptr.numchildren];
		for (int nchild = 0; nchild < treeptr.numchildren; nchild++){
			double[] correc = treeptr.nextptr[nchild].dcorrec;
			NormalDistribution nd = initND(correc);
			for (int nrow = 0; nrow < numrows; nrow++){
				double corr = correc[nrow];
				double pdf = getProportion(nd, corr, minPearson)[1];
				recweight[nrow][nchild] = pdf;
				if(pdf > 0.9){
					IsSeed[nrow] = true;
					if(corr > 0){
						treeptr.nextptr[nchild].posi_genelist.add(nrow);
					}else{
						treeptr.nextptr[nchild].nega_genelist.add(nrow);
					}
				}
			}
		}
		for (int nrow = 0; nrow < numrows; nrow++){
			if(IsSeed[nrow]){
				double sum = 0;
				for (int nchild = 0; nchild < treeptr.numchildren; nchild++){
					sum += recweight[nrow][nchild];
				}
				for (int nchild = 0; nchild < treeptr.numchildren; nchild++){
					recweight[nrow][nchild] /= sum;
				}
			}
		}
		
		if(busego && pathwayData != null && pathwayData.size()>0 && treeptr.numchildren>1){
			double[][][] valuelist = new double[pathwayData.size()][][];
			for(int i=0;i<pathwayData.size();i++) {
	    		RegulatorBindingData knownmodule = pathwayData.get(i);
	    		GettfTable got = new GettfTable(knownmodule, IsSeed, numrows, 400, 0);
	    		got.multithread(knownmodule, IsSeed, numrows, 400);
	    		int[][] coregulateindex = got.coregulateindex;
	    		double[][] coregulatevalue = got.coregulatevalue;
	    		double[][] ptrantable = new double[numrows][treeptr.numchildren];
	    		for (int nrow = 0; nrow < numrows; nrow++){
					if(!IsSeed[nrow]){
						ptrantable[nrow] = FastLogistic2(coregulateindex[nrow], coregulatevalue[nrow], 
								recweight, treeptr.numchildren);
					}
				}
	    		valuelist[i] = ptrantable;
			}
			for (int nrow = 0; nrow < numrows; nrow++){
				if(!IsSeed[nrow]){
					double bestlike = 0;
					double bestcorr = 0;
					int bestchild = -1;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++){
						double corr = treeptr.nextptr[nchild].dcorrec[nrow];
						double pro = getProportion(Normalcorr, corr, minPearson)[0];
						for(int core = 0;core<valuelist.length;core++) {
							pro = pro * valuelist[core][nrow][nchild];
						}
						if(pro > bestlike){
							bestlike = pro;
							bestcorr = corr;
							bestchild = nchild;
						}
					}
					if(bestlike > 0){
						if(bestcorr > 0){
							treeptr.nextptr[bestchild].posi_genelist.add(nrow);
						}else{
							treeptr.nextptr[bestchild].nega_genelist.add(nrow);
						}
					}
				}
			}
		} else {
			for (int nrow = 0; nrow < numrows; nrow++){
				if(!IsSeed[nrow]){
					double bestlike = 0;
					double bestcorr = 0;
					int bestchild = -1;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++){
						double corr = treeptr.nextptr[nchild].dcorrec[nrow];
						double pro = getProportion(Normalcorr, corr, minPearson)[0];
						if(pro > bestlike){
							bestlike = pro;
							bestcorr = corr;
							bestchild = nchild;
						}
					}
					if(bestlike > 0){
						if(bestcorr > 0){
							treeptr.nextptr[bestchild].posi_genelist.add(nrow);
						}else{
							treeptr.nextptr[bestchild].nega_genelist.add(nrow);
						}
					}
				}
			}
		}
	}
	
	public double[] FastLogistic2(int[] coregulatorindex, double[] coregulator,
			double[][] dtrainweight, int numclasses) {
			double ddenom = 0;
			double[] dist = new double[numclasses];
			if(coregulatorindex != null) {
				for (int nclass = 0; nclass < numclasses; nclass++) {
					double dip = 0;
					for (int nindex = 0; nindex < coregulatorindex.length; nindex++){
						int nxindex = coregulatorindex[nindex];
						dip += coregulator[nindex] * dtrainweight[nxindex][nclass];
					}
					double dtempval = Math.exp(dip);
					ddenom += dtempval;
					dist[nclass] = dtempval;
				}
				for(int nindex = 0; nindex < dist.length; nindex++){
					dist[nindex] /= ddenom;
			    }
			}else {
				dist = CONSTANTA(numclasses);
			}
		    return dist;
	}
	public double[] CONSTANTA(int nclass) {
		double[] ptran = new double[nclass];
		for(int i=0;i<ptran.length;i++){
			ptran[i] = (double)1/nclass;
		}
		return ptran;
	}

	public void initAE(Casenode root) {
		if (root != null) {
			for (int nchild = 0; nchild < root.numchildren; nchild++) {
				initAE(root.nextptr[nchild]);
			}
			root.dcorrec = new double[numrows];
			root.dEsum = new double[numcols];
			root.dPsum = new double[numcols];
		}
	}
	
	public void initND(Casenode root) {
		Sumcorr = 0;   
		Sumcorrsq = 0;
		int NDcount = numrows*root.numchildren;
		for (int nrow = 0; nrow < numrows; nrow++){
			 getcorr(traindata[nrow], root, nrow);
		}
		NDmean = Sumcorr / NDcount;
		double dval = (double)NDcount / (double)(NDcount - 1)
				* (Sumcorrsq / NDcount - Math.pow(NDmean, 2));
		NDsigma = Math.sqrt(dval);
		Normalcorr = new NormalDistribution(NDmean, NDsigma);
		//System.out.println(Normalcorr.getMean()+" "+Normalcorr.getStandardDeviation());
	}
	
	public NormalDistribution initND(double[] correc) {	
		Sumcorr = 0;   
		Sumcorrsq = 0;
		int rows = correc.length;
		for (int nrow = 0; nrow < rows; nrow++){
			Sumcorr += correc[nrow];
			Sumcorrsq += Math.pow(correc[nrow], 2);
		}
		NDmean = Sumcorr / rows;
		double dval = (double)rows / (double)(rows - 1)
				* (Sumcorrsq / rows - Math.pow(NDmean, 2));
		NDsigma = Math.sqrt(dval);
		NormalDistribution nd = new NormalDistribution(NDmean, NDsigma);
		return nd;
	}
	
	public void getcorr(double[] vals, Casenode node, int instanceindex){
		Pearsonr pr = new Pearsonr();
		for (int nchild = 0; nchild < node.numchildren; nchild++){
			double[] curve = node.nextptr[nchild].curve;
			if(curve != null){
				double corr = pr.cosineSimilarity(vals, curve);
				Sumcorr += corr;
				Sumcorrsq += Math.pow(corr, 2);
				node.nextptr[nchild].dcorrec[instanceindex] = corr;
			}
		}
	}
	
	private double[] getProportion(NormalDistribution normal, double corr, double traincorr){
		double[] result = new double[2];
		if(Math.abs(corr) > 1){
			throw new IllegalArgumentException("The correlation coefficient should be between -1 and 1.");
		}else if(Math.abs(corr) > traincorr){
			double pdfposi = normal.cumulativeProbability(Math.abs(corr));
			double pdfnega = normal.cumulativeProbability(-Math.abs(corr));
			double pdfdelta = normal.cumulativeProbability(traincorr)
					- normal.cumulativeProbability(-traincorr);
			double P1 = 1 - normal.cumulativeProbability(1);
			double P2 = normal.cumulativeProbability(-1);
			double pdf1 = (pdfposi - pdfnega)/(1-P1-P2);
			double pdf2 = (pdfposi - pdfnega - pdfdelta)/(1-P1-P2-pdfdelta);
			result[0] = pdf1;
			result[1] = pdf2;
			//result = 1 / (1 + Math.exp(-2*Math.log(fdr)));
		}else{
			double pdfposi = normal.cumulativeProbability(Math.abs(corr));
			double pdfnega = normal.cumulativeProbability(-Math.abs(corr));
			double P1 = 1 - normal.cumulativeProbability(1);
			double P2 = normal.cumulativeProbability(-1);
			result[0] = (pdfposi - pdfnega)/(1-P1-P2);
		}
		return result;
	}
	
	public void instanceAE(Casenode root, int bestchild, double[] vals, int nrec) throws Exception{
			double corr = root.nextptr[bestchild].dcorrec[nrec];
			if(corr>0){
				for(int i=0;i<vals.length;i++){
					root.nextptr[bestchild].dEsum[i] += vals[i];
					root.nextptr[bestchild].dPsum[i] += 1;
				}
			}else{
				for(int i=0;i<vals.length;i++){
					root.nextptr[bestchild].dEsum[i] -= vals[i];
					root.nextptr[bestchild].dPsum[i] += 1;
				}
			}
	}
	
	/**
	 * 更新curve
	 */
	public void updateParams(Casenode root) throws Exception {
		for (int nchild = 0; nchild < root.numchildren; nchild++) {
			Casenode ptr = root.nextptr[nchild];
			for(int i=0;i<numcols;i++){
				ptr.curve[i] = ptr.dEsum[i]/ptr.dPsum[i];
			}
		}
	}
	
	/**
	 * 实现baum-welch算法训练模型参数
	 */
	public double trainhmm(Casenode treehmm)throws Exception {
		double dlike = 0;
		double dbestlike = 0;
		initAE(treehmm);
		initND(treehmm);
		
		Gettfbacklg backlg = new Gettfbacklg();
		backlg.multithread(treehmm, Normalcorr, minPearson, numrows, 400);
		int[] bestchildrec = backlg.bestchildrec;
		double[] bestlikerec = backlg.bestlikerec;
		for (int nrow = 0; nrow < numrows; nrow++){
			if(bestlikerec[nrow] > 0){
				instanceAE(treehmm, bestchildrec[nrow], traindata[nrow], nrow);
				dlike += bestlikerec[nrow];
			}
		}
		dbestlike = dlike;
		updateParams(treehmm);
		
		dlike = 0;
		initAE(treehmm);
		initND(treehmm);
		backlg.multithread(treehmm, Normalcorr, minPearson, numrows, 400);
		bestchildrec = backlg.bestchildrec;
		bestlikerec = backlg.bestlikerec;
		for (int nrow = 0; nrow < numrows; nrow++){
			if(bestlikerec[nrow] > 0){
				instanceAE(treehmm, bestchildrec[nrow], traindata[nrow], nrow);
				dlike += bestlikerec[nrow];
			}
		}
		
		int count = 0;
		while((Math.abs(dlike - dbestlike) / dlike) > BEPSILON){
			dbestlike = dlike;
			updateParams(treehmm);
			
			dlike = 0;
			initAE(treehmm);
			initND(treehmm);
			backlg.multithread(treehmm, Normalcorr, minPearson, numrows, 400);
			bestchildrec = backlg.bestchildrec;
			bestlikerec = backlg.bestlikerec;
			for (int nrow = 0; nrow < numrows; nrow++){
				if(bestlikerec[nrow] > 0){
					instanceAE(treehmm, bestchildrec[nrow], traindata[nrow], nrow);
					dlike += bestlikerec[nrow];
				}
			}
			
			count++;
			if(count > 300){
				System.out.println("count>300");
				break;
			}
		}
		
		System.out.println(count+" "+dlike);
		return (dlike);
	}
	
	public void trainhmm2(Casenode treehmm)throws Exception {
		double dlike = 0;
		initAE(treehmm);
		initND(treehmm);
		Gettfbacklg backlg = new Gettfbacklg();
		backlg.multithread(treehmm, Normalcorr, minPearson, numrows, 400);
		//int[] bestchildrec = backlg.bestchildrec;
		double[] bestlikerec = backlg.bestlikerec;
		for (int nrow = 0; nrow < numrows; nrow++){
			if(bestlikerec[nrow] > 0){
				//instanceAE(treehmm, bestchildrec[nrow], traindata[nrow], nrow);
				dlike += bestlikerec[nrow];
			}
		}
		//updateParams(treehmm);
		System.out.println("Final like = " + dlike);
	}
}
