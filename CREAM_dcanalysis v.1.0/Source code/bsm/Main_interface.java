package bsm;


import bsm.core.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import org.apache.commons.math3.distribution.NormalDistribution;

import ssclust.Pearsonr;

import java.util.*;
import java.util.List;
import java.text.*;
import java.net.*;

/**
 * Class implementing the main input interface
 */
public class Main_interface extends JFrame implements ActionListener {
	
	static final String SZDELIM = "|;,";
	static final boolean BDEBUG = false;
	
	BSM_DataSet controldata, casedata;
	RegulatorBindingData pathwayData;
	RegulatorBindingData bindingData;
	RegulatorBindingData mirnaData;
	
	// Main
	static String szoptionalFile = " Optional";
	static String szDataFileDEF = "";
	static boolean bspotcheckDEF = false;
	static String szGeneAnnotationFileDEF = "";
	static String szCrossRefFileDEF = "";
	static int ndbDEF = 0;
	static int nstaticsourceDEF = 0;
	static int numchildDEF = 100;
	static int nummissing = 0;

	// Repeat
	static Vector<String> vRepeatFilesDEF = new Vector<String>();

	// Search Options
	static double dCONVERGENCEDEF = 0.01;
	static double dMinScoreDEF = 0.0;
	static double dDELAYPATHDEF = .15;
	static double dDMERGEPATHDEF = .15;
	static double dPRUNEPATHDEF = .15;
	static double dTHRESHOLD = 0.000001;
	static double dPvalueshold = 0.01;

	// Filtering
	JRadioButton log1, log2, log1_case, log2_case;
	static int nnormalizeDEF = 1;
	static int nnormalizeDEF_case = 1;
	static int ngoDEF = 1;
	static int ntfDEF = 1;
	static int nmirnaDEF = 1;
	ButtonGroup normGroup = new ButtonGroup();
	ButtonGroup normGroup_case = new ButtonGroup();
	
	
	String labels1[] = {" Min value", " Mean value", " Zero"};
	JComboBox comboBox1 = new JComboBox(labels1);
	JComboBox comboBox1_case = new JComboBox(labels1);
	String labels3[] = {" Max value", " Mean value"};
	JComboBox comboBox3 = new JComboBox(labels3);
	JComboBox comboBox3_case = new JComboBox(labels3);
	//String labels4[] = {" Yes"," No"};
	//JComboBox comboBox4 = new JComboBox(labels4);
	JLabel filterlable;
	static int missinter = 0; // Set the missing value as
	static int missinter_case = 0; // Set the missing value as
	static int filterduplicates = 0; // Set the duplicate genes as
	static int filterduplicates_case = 0; // Set the duplicate genes as
	static String szPrefilteredDEF = "";
	static int nMaxMissingDEF = 0;
	static double dMinExpressionDEF = 1;
	static double dMinCorrelationRepeatsDEF = 0;

	// GUI
	static long s1;
	static boolean bendsearch = false;
	static String szcontrolval;
	static String szcaseval;
	static String szmoduleval;
	
	static String missing1, missing_case;
	static boolean btakelog1, btakelog_case;
	

	int npathwaysourcecb;
	int ntfsourcecb;
	int nppisourcecb;
	String szuserFileField1;
	String szuserFileField2;
	String szuserFileField3;
	
	static JFileChooser theChooser = new JFileChooser();
	
	// Regulator Scoring GUI
	JButton regScoreHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JTextField regScoreField;
	static int NUMCOLS = 42;
	// Strings for the labels
	static Color gray = new Color(235,235,235);

	JButton infoButton = new JButton(Util.createImageIcon("About16.gif"));
	static JFileChooser fc1 = new JFileChooser(new File("Enrichment_analysis_files"));
	static JFileChooser fc2 = new JFileChooser(new File("Expression data (sample)"));
	static JFileChooser fc3 = new JFileChooser(new File("Module files (sample)"));

	static Container contentPane;
	JPanel textpanel;
	JPanel panel4;

	JButton controlbutton= new JButton();
	JButton casebutton= new JButton();
	JButton modulebutton= new JButton();
	JTextArea runningtext;
	JButton Run = new JButton();
	JButton currentButton = new JButton();
	JTextField controlField, caseField, moduleField;
	final JSpinner j18, j18_case;
	HashMap<String, Integer> geneindex = new HashMap<String, Integer>();
	HashMap<String, List<Integer>> resultlist, posilist, negalist;

	/**
	 * Class constructor - builds the input interface calls parseDefaults to get
	 * the initial settings from a default settings file if specified
	 */
	public Main_interface() throws FileNotFoundException, IOException {
		super("CREAM_DCanalysis v.1.0");

		contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
		contentPane.setLayout(layout);

		JPanel panel1=new JPanel();
		panel1.setBorder(new TitledBorder(null,"Load Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel1.setLayout(new BorderLayout());
		
		JLabel j7=new JLabel("Load expression data 1:");
		j7.setBounds(15,20,300,35);
		controlField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		controlField.setBounds(170,23,400,25);
		controlField.setBorder(BorderFactory.createLineBorder(Color.black));
		controlField.setOpaque(true);
		controlField.setBackground(Color.white);
		controlbutton.setText("Load File");
		controlbutton.setHideActionText(true);
		controlbutton.addActionListener(this);
		controlbutton.setBounds(575,23,90,25);
		
		/*****************************************************************/
		
		JLabel j8=new JLabel("Load expression data 2:");
		j8.setBounds(15,55,300,35);
		caseField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		caseField.setBounds(170,58,400,25);
		caseField.setBorder(BorderFactory.createLineBorder(Color.black));
		caseField.setOpaque(true);
		caseField.setBackground(Color.white);
		casebutton.setText("Load File");
		casebutton.setHideActionText(true);
		casebutton.addActionListener(this);
		casebutton.setBounds(575,58,90,25);
		
		/*****************************************************************/
		
		JLabel j9=new JLabel("Load module file:");
		j9.setBounds(15,90,300,35);
		moduleField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		moduleField.setBounds(170,93,400,25);
		moduleField.setBorder(BorderFactory.createLineBorder(Color.black));
		moduleField.setOpaque(true);
		moduleField.setBackground(Color.white);
		modulebutton.setText("Load File");
		modulebutton.setHideActionText(true);
		modulebutton.addActionListener(this);
		modulebutton.setBounds(575,93,90,25);
		
		/*****************************************************************/

		panel1.setLayout(null);
		panel1.add(j7);
		panel1.add(j8);
		panel1.add(j9);

		panel1.add(controlField);
		panel1.add(caseField);
		panel1.add(moduleField);
		
		panel1.add(controlbutton);
		panel1.add(casebutton);
		panel1.add(modulebutton);
		
        JScrollPane   scrollpanel1   =   new   JScrollPane(panel1, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel1.setBounds(100, 100, 745, 70);
		panel1.setPreferredSize(new Dimension(scrollpanel1.getWidth() - 50, scrollpanel1.getHeight()*2));
		
		/*******************************************************/
		
		JPanel panel3=new JPanel();
		JScrollPane   scrollpanel2   =   new   JScrollPane(panel3, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel2.setBounds(100, 100, 745, 100);
        panel3.setPreferredSize(new Dimension(scrollpanel2.getWidth() - 50, scrollpanel2.getHeight()*2));
		panel3.setBorder(new TitledBorder(null,"Data Preprocessing",TitledBorder.LEFT,TitledBorder.TOP));
		panel3.setLayout(new BorderLayout());
		
		JPanel panel31=new JPanel();
		panel31.setBorder(new TitledBorder(null,"Control Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel31.setBounds(5,20,365,175);
		panel31.setLayout(null);
		
		JLabel log=new JLabel("Log normalize data:");
		log.setBounds(15,20,300,35);
	    log1 = new JRadioButton("Yes");
		log2 = new JRadioButton("No");
		log1.setBounds(140, 25, 50, 22);
		log2.setBounds(190, 25, 50, 22);
		if (nnormalizeDEF == 0) {
			log1.setSelected(true);
		} else if (nnormalizeDEF == 1) {
			log2.setSelected(true);
		}
		
		JLabel mm=new JLabel("Max number of missing values:");
		mm.setBounds(15,55,300,35);
		SpinnerNumberModel misscontrol = new SpinnerNumberModel(new Integer(nummissing), new Integer(0), null, new Integer(1));
		j18 = new JSpinner(misscontrol);
		j18.setBounds(210,60,45,22);
		
		JLabel inter=new JLabel("Set missing values as:");
		inter.setBounds(15,95,300,35);
		comboBox1.setBounds(155,100,100,21);
	    if (missinter == 0) {
	    	comboBox1.setSelectedIndex(0);
		} else if(missinter == 1){
			comboBox1.setSelectedIndex(1);
		} else{
			comboBox1.setSelectedIndex(2);
		}
	    
	    JLabel duplicates=new JLabel("Set duplicate genes as:");
		duplicates.setBounds(15,135,300,35);
	    comboBox3.setBounds(155,140,100,21);
	    if (filterduplicates == 0) {
	    	comboBox3.setSelectedIndex(0);
		} else if(filterduplicates == 1){
			comboBox3.setSelectedIndex(1);
		}
		
		panel31.add(log);
		panel31.add(log1);
		panel31.add(log2);
		normGroup.add(log1);
		normGroup.add(log2);
		panel31.add(mm);
		panel31.add(j18);
		panel31.add(inter);
		panel31.add(comboBox1);
		panel31.add(duplicates);
		panel31.add(comboBox3);
		
		/*******************************************************/
	    JPanel panel32=new JPanel();
		panel32.setBorder(new TitledBorder(null,"Case Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel32.setBounds(370,20,365,175);
		panel32.setLayout(null);

		JLabel log_case=new JLabel("Log normalize data:");
		log_case.setBounds(15,20,300,35);
	    log1_case = new JRadioButton("Yes");
		log2_case = new JRadioButton("No");
		log1_case.setBounds(140, 25, 50, 22);
		log2_case.setBounds(190, 25, 50, 22);
		if (nnormalizeDEF_case == 0) {
			log1_case.setSelected(true);
		} else if (nnormalizeDEF_case == 1) {
			log2_case.setSelected(true);
		}
		
		JLabel mm_case=new JLabel("Max number of missing values:");
		mm_case.setBounds(15,55,300,35);
		SpinnerNumberModel misscontrol_case = new SpinnerNumberModel(new Integer(nummissing), new Integer(0), null, new Integer(1));
		j18_case = new JSpinner(misscontrol_case);
		j18_case.setBounds(210,60,45,22);
		
		JLabel inter_case=new JLabel("Set missing values as:");
		inter_case.setBounds(15,95,300,35);
		comboBox1_case.setBounds(155,100,100,21);
	    if (missinter_case == 0) {
	    	comboBox1_case.setSelectedIndex(0);
		} else if(missinter_case == 1){
			comboBox1_case.setSelectedIndex(1);
		} else{
			comboBox1_case.setSelectedIndex(2);
		}
	    
	    JLabel duplicates_case=new JLabel("Set duplicate genes as:");
		duplicates_case.setBounds(15,135,300,35);
	    comboBox3_case.setBounds(155,140,100,21);
	    if (filterduplicates_case == 0) {
	    	comboBox3_case.setSelectedIndex(0);
		} else if(filterduplicates_case == 1){
			comboBox3_case.setSelectedIndex(1);
		}

	    panel32.add(log_case);
		panel32.add(log1_case);
		panel32.add(log2_case);
		normGroup_case.add(log1_case);
		normGroup_case.add(log2_case);
		panel32.add(mm_case);
		panel32.add(j18_case);
		panel32.add(inter_case);
		panel32.add(comboBox1_case);
		panel32.add(duplicates_case);
		panel32.add(comboBox3_case);
		
		panel3.setLayout(null);
		panel3.add(panel31); 
		panel3.add(panel32);
		/********************************************************/
		
		textpanel=new JPanel();
		textpanel.setPreferredSize(new Dimension(745, 150));
		textpanel.setBorder(new TitledBorder(null,"Analysis",TitledBorder.LEFT,TitledBorder.TOP));
		BoxLayout layout2 = new BoxLayout(textpanel, BoxLayout.X_AXIS);
		textpanel.setLayout(layout2);
		
		panel4=new JPanel();
		panel4.setLayout(new BorderLayout());
		runningtext=new JTextArea();
		runningtext.setLineWrap(true);
		JScrollPane sp1=new JScrollPane(runningtext);
		panel4.add(sp1);
		textpanel.add(panel4);

		Run.setText("Run");
		Run.setHideActionText(true);
		Run.addActionListener(this);
		infoButton.addActionListener(this);
		infoButton.setBackground(gray);
		
		JPanel panel6=new JPanel();
		panel6.setPreferredSize(new Dimension(745, 50));
		panel6.add(Run);
		panel6.add(infoButton);
		contentPane.add(scrollpanel1);
		contentPane.add(scrollpanel2);
		contentPane.add(textpanel);
		contentPane.add(panel6);
	}

	/**
	 * Checks if the two data sets have the same number of rows, time points,
	 * and the gene name matches.
	 */
	public static void errorcheck(BSM_DataSet theDataSet1,
			BSM_DataSet theOtherSet) {

		if (theDataSet1.numcols != theOtherSet.numcols) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of columns as original, expecting "
							+ theDataSet1.numcols + " found "
							+ theOtherSet.numcols + " in the repeat");
		} else if (theDataSet1.numrows != theOtherSet.numrows) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of spots as the original, expecting "
							+ theDataSet1.numrows + " found "
							+ theOtherSet.numrows + " in the repeat");
		} else {
			for (int nrow = 0; nrow < theDataSet1.numrows; nrow++) {
				if (!theDataSet1.genenames[nrow]
						.equals(theOtherSet.genenames[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ theDataSet1.genenames[nrow] + " found "
							+ theOtherSet.genenames[nrow]);
				} else if (!theDataSet1.probenames[nrow]
						.equals(theOtherSet.probenames[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ theDataSet1.probenames[nrow] + " found "
							+ theOtherSet.probenames[nrow]);
				}
			}
		}
	}

	/**
	 * Checks if origcols and nrepeat cols are the same value, the length of
	 * origgenes and repeatgenes is the same, and the gene names are the same
	 */
	public static void errorcheck(String[] origgenes, String[] repeatgenes,
			int norigcols, int nrepeatcols) {
		if (norigcols != nrepeatcols) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of columns as original, expecting "
							+ norigcols + " found " + nrepeatcols
							+ " in the repeat");
		} else if (origgenes.length != repeatgenes.length) {
			throw new IllegalArgumentException(
					"Repeat data set must have same "
							+ "number of spots as the original, expecting "
							+ origgenes.length + " found " + repeatgenes.length
							+ " in the repeat");
		} else {
			for (int nrow = 0; nrow < origgenes.length; nrow++) {
				if (!origgenes[nrow].equals(repeatgenes[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ origgenes[nrow] + " found " + repeatgenes[nrow]);
				} else if (!origgenes[nrow].equals(repeatgenes[nrow])) {
					throw new IllegalArgumentException("In row " + nrow
							+ " of the repeat set " + "expecting gene symbol "
							+ origgenes[nrow] + " found " + repeatgenes[nrow]);
				}
			}
		}
	}

	private BSM_DataSet buildset(String szexp1val, 
			Integer maxcase, Integer missinter, boolean btakelog, int filterduplicates) 
			throws Exception {
		BSM_DataSet theDataSet1 = new BSM_DataSet(szexp1val, maxcase, missinter, btakelog);
		theDataSet1 = new BSM_DataSet(theDataSet1.filterMissing());	
		if(filterduplicates==0){
			theDataSet1 = new BSM_DataSet(theDataSet1.maxAndFilterDuplicates());
		}else if(filterduplicates==1){
			theDataSet1 = new BSM_DataSet(theDataSet1.averageAndFilterDuplicates());
		}
		return theDataSet1;	
	}
	
	/**
	 * A control method that handles the response for when the execute button on
	 * the interface is pressed including building the data set, running the
	 * TSMiner modeling procedure, and displaying the results
	 */
	public void clusterscript(String szcontrolval, String szcaseval,
			String szmoduleval, String missing1, int missinter1, boolean btakelog1, int filterduplicates,
			String missing2, int missinter2, boolean btakelog2, int filterduplicates2) throws Exception {
		
		//szcontrolval = "D:\\Personal issues\\CREAM\\CREAM tool\\CREAM_dcanalysis v.1.0\\Expression data (sample)/LIHC_stage4.txt";
		//szcaseval = "D:\\Personal issues\\CREAM\\CREAM tool\\CREAM_dcanalysis v.1.0\\Expression data (sample)/LIHC_stage2.txt";
		//szmoduleval = "D:\\Personal issues\\CREAM\\CREAM tool\\CREAM_dcanalysis v.1.0\\Module files (sample)/LIHCstage4.txt";
		
		List<File> go = getAllFile(new File("Enrichment_analysis_files"));
		
				if (szcontrolval.trim().equals("") ) {
					throw new IllegalArgumentException("No caontrol expression data given!");
				} else if (!(new File(szcontrolval)).exists()) {
					throw new IllegalArgumentException("The control data file '" + szcontrolval+ "' cannot be found.");
				} else if (szcaseval.trim().equals("")) {
					throw new IllegalArgumentException("No case expression data given!");
				} else if (!(new File(szcaseval)).exists()) {
					throw new IllegalArgumentException("The case data file '" + szcaseval+ "' cannot be found.");
				}else if (szmoduleval.trim().equals("")) {
					throw new IllegalArgumentException("No module file given!");
				} else if (!(new File(szmoduleval)).exists()) {
					throw new IllegalArgumentException("The module file '" + szcaseval+ "' cannot be found.");
				}
				
				controldata = buildset(szcontrolval, Integer.parseInt(missing1), missinter1, btakelog1, filterduplicates);
				controldata.normalization();
				casedata = buildset(szcaseval, Integer.parseInt(missing2), missinter2, btakelog2, filterduplicates2);
				casedata.normalization();
				
				for(int i=0;i<controldata.genenames.length;i++){
					geneindex.put(controldata.genenames[i].toUpperCase(), i);
				}
				System.out.println(geneindex.size());
				
				if(go != null && go.size()>0){
					pathwayData = new RegulatorBindingData(go.get(0).toString(),
							controldata, SZDELIM, null, null);
					pathwayData.deDuplication();
				}
				
				readModuleFile(szmoduleval);
				
				try {
					DCanalysis caseModules = new DCanalysis(controldata, casedata, pathwayData, 
							resultlist, posilist, negalist, runningtext);
		        	final DCanalysis finalcase = caseModules;
		        	
		        	javax.swing.SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							JFrame frame1 = new JFrame("DC Analysis Results");
							Container theDialogContainer = frame1.getContentPane();
							theDialogContainer.setBackground(Color.white);
							JTabbedPane tabbedPane = new JTabbedPane();
								CaseTable newContentPane1 = new CaseTable(frame1, 
										finalcase, finalcase.ModuleSet, controldata, casedata);
								newContentPane1.setOpaque(true); // content panes must be opaque
								tabbedPane.addTab("Resulting Modules", null, newContentPane1,"TF Enrichment");
							
							theDialogContainer.add(tabbedPane);
							frame1.setPreferredSize(new Dimension(800, 800));
							frame1.setMinimumSize(new Dimension(800, 800));
							frame1.setContentPane(theDialogContainer);
							frame1.pack();
							frame1.setVisible(true);
						}
					});
					
				} catch (Exception e) {
					e.printStackTrace();
			}
		
		//System.exit(0);
		long e1 = System.currentTimeMillis();
		System.out.println("Time: " + (e1 - s1) + "ms");
	}
	
	public List<File> getAllFile(File dirFile) {
        // 如果文件夹不存在或着不是文件夹，则返回 null
        if (Objects.isNull(dirFile) || !dirFile.exists() || dirFile.isFile())
            return null;

        File[] childrenFiles = dirFile.listFiles();
        if (Objects.isNull(childrenFiles) || childrenFiles.length == 0)
            return null;

        List<File> files = new ArrayList<>();
        for (File childFile : childrenFiles) {
            if (childFile.isFile()) {
                files.add(childFile);
            }
        }
        return files;
    }
	
	public List<List<String>> readFile(String filename){
        List<List<String>> value = new ArrayList<List<String>>();
		try {
			BufferedReader bufferedReader1 = new BufferedReader(new FileReader(filename));   
			String lineTxt = null;
            while((lineTxt = bufferedReader1.readLine()) != null){
                value.add(Arrays.asList(lineTxt.split("\t")));
            }
            bufferedReader1.close();  
                
    } catch (Exception e) {
    	throw new IllegalArgumentException("Read saved model file error.");
    }
		return value;
	}
	
	public void readModuleFile(String filename){
		resultlist = new HashMap<String, List<Integer>>();
		posilist = new HashMap<String, List<Integer>>();
		negalist = new HashMap<String, List<Integer>>();
		try {
			BufferedReader bufferedReader1 = new BufferedReader(new FileReader(filename));   
			String lineTxt = bufferedReader1.readLine();
			while((lineTxt = bufferedReader1.readLine()) != null){
            	if(Arrays.asList(lineTxt.split("\t")).size() > 2){
            		String moduleid = Arrays.asList(lineTxt.split("\t")).get(0);
            		String genename = Arrays.asList(lineTxt.split("\t")).get(1).toUpperCase(Locale.ENGLISH);
            		int correlation = Integer.parseInt(Arrays.asList(lineTxt.split("\t")).get(2));
            		if(geneindex.get(genename) != null) {
            			if(correlation > 0) {
            				if(posilist.get(moduleid) != null) {
                				List<Integer> list = posilist.get(moduleid);
                				if(!list.contains(geneindex.get(genename))) {
                					list.add(geneindex.get(genename));
                				}
                				posilist.put(moduleid, list);
                    		} else {
                        		List<Integer> list = new ArrayList<Integer>();
                        		list.add(geneindex.get(genename));
                        		posilist.put(moduleid, list);
                    		}
            			} else if(correlation < 0) {
            				if(negalist.get(moduleid) != null) {
                				List<Integer> list = negalist.get(moduleid);
                				if(!list.contains(geneindex.get(genename))) {
                					list.add(geneindex.get(genename));
                				}
                				negalist.put(moduleid, list);
                    		} else {
                        		List<Integer> list = new ArrayList<Integer>();
                        		list.add(geneindex.get(genename));
                        		negalist.put(moduleid, list);
                    		}
            			}
            			if(resultlist.get(moduleid) != null) {
            				List<Integer> list = resultlist.get(moduleid);
            				if(!list.contains(geneindex.get(genename))) {
            					list.add(geneindex.get(genename));
            				}
            				resultlist.put(moduleid, list);
                		} else {
                    		List<Integer> list = new ArrayList<Integer>();
                    		list.add(geneindex.get(genename));
                    		resultlist.put(moduleid, list);
                		}
            			
            		}
            	} else if(Arrays.asList(lineTxt.split("\t")).size() == 2){
            		String moduleid = Arrays.asList(lineTxt.split("\t")).get(0);
            		String genename = Arrays.asList(lineTxt.split("\t")).get(1).toUpperCase(Locale.ENGLISH);
            		if(geneindex.get(genename) != null) {
            			if(posilist.get(moduleid) != null) {
                			List<Integer> list = posilist.get(moduleid);
                			if(!list.contains(geneindex.get(genename))) {
                				list.add(geneindex.get(genename));
                			}
                			posilist.put(moduleid, list);
                    	} else {
                        	List<Integer> list = new ArrayList<Integer>();
                        	list.add(geneindex.get(genename));
                        	posilist.put(moduleid, list);
                    	}
            			if(resultlist.get(moduleid) != null) {
            				List<Integer> list = resultlist.get(moduleid);
            				if(!list.contains(geneindex.get(genename))) {
            					list.add(geneindex.get(genename));
            				}
            				resultlist.put(moduleid, list);
                		} else {
                    		List<Integer> list = new ArrayList<Integer>();
                    		list.add(geneindex.get(genename));
                    		resultlist.put(moduleid, list);
                		}
            			
            		}
            	} else {
            		System.out.println(lineTxt);
            	}
            }
            bufferedReader1.close();
    } catch (Exception e) {
    	throw new IllegalArgumentException("Read module file error.");
    }	
	}

	/**
	 * define the button methods
	 */
	public void actionPerformed(ActionEvent e) {
		Object esource = e.getSource();

		if (esource == controlbutton) {
			int returnVal = fc2.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc2.getSelectedFile();
				controlField.setText(file.getAbsolutePath());
			}
		} else if (esource == casebutton) {
			int returnVal = fc2.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc2.getSelectedFile();
				caseField.setText(file.getAbsolutePath());
			}
		} else if (esource == modulebutton) {
			int returnVal = fc3.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc3.getSelectedFile();
				moduleField.setText(file.getAbsolutePath());
			}
		} else if (esource == Run) {
			s1 = System.currentTimeMillis();
			szcontrolval = controlField.getText();
			szcaseval = caseField.getText(); 
			szmoduleval = moduleField.getText(); 
			
			missing1 = j18.getValue().toString();
			missing_case = j18_case.getValue().toString();
			
			int temp2=comboBox1.getSelectedIndex();
			if(temp2==0){
				missinter=0;
			}else if(temp2==1){
				missinter=1;
			}else{
				missinter=2;
			}
			temp2=comboBox1_case.getSelectedIndex();
			if(temp2==0){
				missinter_case=0;
			}else if(temp2==1){
				missinter_case=1;
			}else{
				missinter_case=2;
			}
			
			int temp3=comboBox3.getSelectedIndex();
			if(temp3==0){
				filterduplicates=0;
			}else if(temp3==1){
				filterduplicates=1;
			}
			temp3=comboBox3_case.getSelectedIndex();
			if(temp3==0){
				filterduplicates_case=0;
			}else if(temp3==1){
				filterduplicates_case=1;
			}
			
			btakelog1 = log1.isSelected();
			btakelog_case = log1_case.isSelected();

			this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			final JFrame fframe = this;

			Runnable clusterrun = new Runnable() {
				public void run() {
					Run.setEnabled(false);
					try {
						clusterscript(szcontrolval, szcaseval, szmoduleval, 
								missing1, missinter, btakelog1, filterduplicates, 
								missing_case, missinter_case, btakelog_case, filterduplicates_case);
					} catch (IllegalArgumentException iex) {
						final IllegalArgumentException fiex = iex;
						iex.printStackTrace(System.out);

						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, fiex.getMessage(), "Error",JOptionPane.ERROR_MESSAGE);
							}
						});
					} catch (Exception ex) {
						final Exception fex = ex;

						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, fex.toString(), "Exception thrown",JOptionPane.ERROR_MESSAGE);
								fex.printStackTrace(System.out);
							}
						});
					}
					Run.setEnabled(true); 
					fframe.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
				}
			};
			(new Thread(clusterrun)).start();
			
		} else if (esource == infoButton) {//The Information and Help button
			String szMessage = "This is version 1.0.0 of CREAM.\n\n"
					+ "The CREAM is available under a GPL v3.0 license.\n"
					+ "Any questions or bugs found should "
					+ "be emailed to free1234hm@163.com.";

			Util.renderDialog(this, szMessage, 50, 100, "Information");
		}
	}
	
	/**
	 * Create the GUI and show it. For thread safety, this method should be
	 * invoked from the event-dispatching thread.
	 */
	private static void createAndShowGUI() throws FileNotFoundException,IOException {
		// Make sure we have nice window decorations.
		// JFrame.setDefaultLookAndFeelDecorated(true);
		// Create and set up the window.
		JFrame frame = new Main_interface();
		frame.setLocation(10, 25);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

		// Display the window.
		frame.pack();
		frame.setVisible(true);
	}

	
	/**
	 * The main method which when executed will have the input interface created
	 */
	public static void main(String[] args) throws Exception {
			javax.swing.SwingUtilities.invokeLater(new Runnable() {
				public void run() {
					try {
						createAndShowGUI();
					} catch (FileNotFoundException ex) {
						ex.printStackTrace(System.out);
					} catch (IOException ex) {
						ex.printStackTrace(System.out);
					}
				}
			});
	}
}
