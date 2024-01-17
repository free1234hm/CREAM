package bsm;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;
import bsm.core.DataSetCore;

public class PPIInteractionData {

	HashMap<Integer, List<Integer>> ppiindex;
	HashMap<Integer, List<Double>> ppivalue;
	String REGTYPE = "Protein";
	
	public PPIInteractionData(String tfBindingDataFile,
			BSM_DataSet theds, boolean[]keepgene, int nkeep) throws IOException {

		ppiindex = new HashMap<Integer, List<Integer>>();
		ppivalue = new HashMap<Integer, List<Double>>();
		String dataFiles = tfBindingDataFile;

			if (dataFiles != null && !dataFiles.equals("")) {
				BufferedReader br = null;
				try {//¶ÁÎÄ¼þ
					br = new BufferedReader(new InputStreamReader(
							new GZIPInputStream(new FileInputStream(dataFiles))));
				} catch (IOException ex) {
					br = new BufferedReader(new FileReader(dataFiles));
				}

				String szLine = br.readLine();
				StringTokenizer st = new StringTokenizer(szLine, "\t");
				String szh1 = "";

				int numType;
				if (szLine == null) {
					throw new IllegalArgumentException("Empty PPI interaction input file found!");
				} else if (szLine.startsWith("\t")) {
					numType = st.countTokens();
				} else {
					numType = st.countTokens() - 1;
					szh1 = st.nextToken();
				}
				// Read the top row of the file in      //¶Á¶¥ÐÐ
				String[] tempRegNames = new String[numType];
				for (int nRegIndex = 0; st.hasMoreTokens(); nRegIndex++) {
					tempRegNames[nRegIndex] = st.nextToken();
				}
				
				if(!szh1.equalsIgnoreCase(REGTYPE)){
					throw new IllegalArgumentException("PPI file should start with "
							+ "'Protein' in the first row");
				}

				while ((szLine = br.readLine()) != null) {
					st = new StringTokenizer(szLine, "\t");
					String szprotein1 = st.nextToken().toUpperCase(Locale.ENGLISH);
					String szprotein2 = st.nextToken().toUpperCase(Locale.ENGLISH);
					double ninput;
					if (st.hasMoreTokens()) {
						String szToken = st.nextToken();
						try {
							ninput = Double.parseDouble(szToken);
						} catch (NumberFormatException nfex) {
							throw new IllegalArgumentException(szToken + " is not a"
									+ " valid score for a protein-protein interaction");
						}
					} else {
						ninput = 1.0;
					}
					
					Integer protein1 = theds.gene2int.get(szprotein1);
					Integer protein2 = theds.gene2int.get(szprotein2);
					
					if(protein1 != null && protein2 != null){
						List<Integer> prolist1 = ppiindex.get(protein1);
						List<Double> vallist1 = ppivalue.get(protein1);
						List<Integer> prolist2 = ppiindex.get(protein2);
						List<Double> vallist2 = ppivalue.get(protein2);
							
							if(prolist1 == null){
								prolist1 = new ArrayList<Integer>();
								vallist1 = new ArrayList<Double>();
								prolist1.add(protein2);
								vallist1.add(ninput);
								ppiindex.put(protein1, prolist1);
								ppivalue.put(protein1, vallist1);
							}else{
								if(!prolist1.contains(protein2)){
									prolist1.add(protein2);
									vallist1.add(ninput);
									ppiindex.put(protein1, prolist1);
									ppivalue.put(protein1, vallist1);
								}
							}
							if(prolist2 == null){
								prolist2 = new ArrayList<Integer>();
								vallist2 = new ArrayList<Double>();
								prolist2.add(protein1);
								vallist2.add(ninput);
								ppiindex.put(protein2, prolist2);
								ppivalue.put(protein2, vallist2);
							}else{
								if(!prolist2.contains(protein1)){
									prolist2.add(protein1);
									vallist2.add(ninput);
									ppiindex.put(protein2, prolist2);
									ppivalue.put(protein2, vallist2);
								}
							}
					}
				}
			}
	}
	public int contain(String tfname, BSM_DataSet allds) {
		int index = -1;
		String[] list = allds.genenames;
		for(int i=0;i<list.length;i++){
			if(tfname.equalsIgnoreCase(list[i])){
				index = i;
				break;
			}
		}
		return index;
	}
}
