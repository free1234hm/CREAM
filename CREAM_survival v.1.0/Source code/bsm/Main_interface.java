package bsm;


import bsm.core.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;
import javax.swing.event.*;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.jfree.data.category.DefaultCategoryDataset;

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
	
	BSM_DataSet expressiondata;
	HashMap<String, List<Integer>> posimodules, negamodules;
	RegulatorBindingData pathwayData;
	RegulatorBindingData bindingData;
	RegulatorBindingData mirnaData;
	
	// Main
	static String szoptionalFile = " Optional";
	static String szDataFileDEF = "";
	static int nstaticsourceDEF = 0;
	static int nummissing = 0;

	// GUI
	static long s1;
	static boolean bendsearch = false;
	static String szsurvival;
	static String szdataval;
	static String szmoduleval;
	

	int npathwaysourcecb;
	int ntfsourcecb;
	int nppisourcecb;
	String szuserFileField1;
	String szuserFileField2;
	String szuserFileField3;
	
	static JFileChooser theChooser = new JFileChooser();

	JComboBox survivalsourcecb;
	JComboBox tfsourcecb;
	JComboBox mirnasourcecb;

	// Regulator Scoring GUI
	JButton regScoreFileButton = new JButton("Browse...", Util.createImageIcon("Open16.gif"));
	JButton regScoreHButton = new JButton(Util.createImageIcon("Help16.gif"));
	JTextField regScoreField;
	static int NUMCOLS = 42;
	// Strings for the labels
	static Color gray = new Color(235,235,235);
	static Color defaultColor;
	static String[] staticsourceArray1 = { "User provided" };

	JTextField goField;
	JTextField extraField;
	JTextField xrefField;
	JTextField categoryIDField;
	JTextField taxonField;
	JTextField evidenceField;
	JButton infoButton = new JButton(Util.createImageIcon("About16.gif"));
	static String fc = "TCGA survival data";
	static JFileChooser fc1 = new JFileChooser(new File("TCGA survival data"));
	static JFileChooser fc2 = new JFileChooser(new File("TCGA expression data"));
	static JFileChooser fc3 = new JFileChooser(new File("Module files (sample)"));

	static Container contentPane;
	JPanel textpanel;
	JPanel panel4;


	JButton survivalbutton= new JButton();
	JButton databutton= new JButton();
	JButton modulebutton= new JButton();
	JTextArea runningtext;
	JButton Run = new JButton();
	JTextField survivalField, dataField, moduleField;
	static String missing1;
	static boolean btakelog1;
	static int missinter = 0; // Set the missing value as
	static int filterduplicates = 0; // Set the duplicate genes as
	static String division;
	
	JRadioButton log1, log2;
	static int nnormalizeDEF = 1;
	final JSpinner j18;
	String labels1[] = {" Min value", " Mean value", " Zero"};
	JComboBox comboBox1 = new JComboBox(labels1);
	JComboBox comboBox1_case = new JComboBox(labels1);
	String labels3[] = {" Max value", " Mean value"};
	JComboBox comboBox3 = new JComboBox(labels3);
	JComboBox comboBox3_case = new JComboBox(labels3);
	ButtonGroup normGroup = new ButtonGroup();
	final JSpinner jpro;
	HashMap<String, List<Integer>> resultlist, posilist, negalist;
	HashMap<String, Integer> geneindex;
	

	/**
	 * Class constructor - builds the input interface calls parseDefaults to get
	 * the initial settings from a default settings file if specified
	 */
	public Main_interface() throws FileNotFoundException, IOException {
		super("CREAM_survival v.1.0");

		File dir1 = new File(fc);
		{
			String[] children = dir1.list();
			if (children == null) {
			} else {
				staticsourceArray1 = new String[children.length + 1];
				staticsourceArray1[0] = "User Provided";
				for (int i = 0; i < children.length; i++) {
					// Get filename of file or directory
					staticsourceArray1[i + 1] = children[i];
				}
			}
		}
		survivalsourcecb = new JComboBox(staticsourceArray1);
		survivalsourcecb.addActionListener(this);

		contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
		contentPane.setLayout(layout);

		JPanel panel1=new JPanel();
		panel1.setBorder(new TitledBorder(null,"Load Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel1.setLayout(new BorderLayout());
		
     /*****************************************************************/
		
		JLabel j1=new JLabel("Survival file source:");
		j1.setBounds(15,20,200,35);
		survivalsourcecb.setBounds(170, 25, 150, 25);
		
		JLabel j7=new JLabel("Load survival file:");
		j7.setBounds(15,55,300,35);
		survivalField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		survivalField.setBounds(170,58,400,25);
		survivalField.setBorder(BorderFactory.createLineBorder(Color.black));
		survivalField.setOpaque(true);
		survivalField.setBackground(Color.white);
		survivalbutton.setText("Load File");
		survivalbutton.setHideActionText(true);
		survivalbutton.addActionListener(this);
		survivalbutton.setBounds(575,58,90,25);
		
		/*****************************************************************/
		
		JLabel j8=new JLabel("Load expression data:");
		j8.setBounds(15,90,300,35);
		dataField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		dataField.setBounds(170,93,400,25);
		dataField.setBorder(BorderFactory.createLineBorder(Color.black));
		dataField.setOpaque(true);
		dataField.setBackground(Color.white);
		databutton.setText("Load File");
		databutton.setHideActionText(true);
		databutton.addActionListener(this);
		databutton.setBounds(575,93,90,25);
		
		/*****************************************************************/
		
		JLabel j9=new JLabel("Load module file:");
		j9.setBounds(15,125,300,35);
		moduleField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		moduleField.setBounds(170,128,400,25);
		moduleField.setBorder(BorderFactory.createLineBorder(Color.black));
		moduleField.setOpaque(true);
		moduleField.setBackground(Color.white);
		modulebutton.setText("Load File");
		modulebutton.setHideActionText(true);
		modulebutton.addActionListener(this);
		modulebutton.setBounds(575,128,90,25);
		
		/*****************************************************************/
		
		panel1.setLayout(null);
		panel1.add(j1);
		panel1.add(survivalsourcecb);
		survivalsourcecb.setSelectedIndex(nstaticsourceDEF);

		panel1.add(j7);
		panel1.add(j8);
		panel1.add(j9);
		
		panel1.add(survivalField);
		panel1.add(dataField);
		panel1.add(moduleField);

		panel1.add(survivalbutton);
		panel1.add(databutton);
		panel1.add(modulebutton);
		
        JScrollPane   scrollpanel1   =   new   JScrollPane(panel1, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel1.setBounds(100, 100, 745, 90);
		panel1.setPreferredSize(new Dimension(scrollpanel1.getWidth() - 50, scrollpanel1.getHeight()*2));
		
		/**********************************************************************/
		JPanel panel3=new JPanel();
		JScrollPane   scrollpanel2   =   new   JScrollPane(panel3, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel2.setBounds(100, 100, 745, 100);
        panel3.setPreferredSize(new Dimension(scrollpanel2.getWidth() - 50, scrollpanel2.getHeight()*2));
		panel3.setBorder(new TitledBorder(null,"Set parameters",TitledBorder.LEFT,TitledBorder.TOP));
		panel3.setLayout(new BorderLayout());
		
		JPanel panel31=new JPanel();
		panel31.setBorder(new TitledBorder(null,"Data preprocessing",TitledBorder.LEFT,TitledBorder.TOP));
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

		/******************************************************************/
		JPanel panel32=new JPanel();
		panel32.setBorder(new TitledBorder(null,"Survival analysis",TitledBorder.LEFT,TitledBorder.TOP));
		panel32.setBounds(370,20,365,175);
		panel32.setLayout(null);
		
		JLabel ll1=new JLabel("Compare the top and bottom");
		ll1.setBounds(15,20,300,35);
		SpinnerNumberModel firstgroup = new SpinnerNumberModel(0.25, 0, 0.5, 0.01);
		jpro = new JSpinner(firstgroup);
		jpro.setBounds(185,25,45,22);
		JLabel ll2=new JLabel("of patients.");
		ll2.setBounds(235,20,300,35);
		
		panel32.add(ll1);
		panel32.add(jpro);
		panel32.add(ll2);
		
		panel3.setLayout(null);
		panel3.add(panel31);
		panel3.add(panel32);
		/***************************************************************/
		
		
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
	public void clusterscript(String szsurvival, String szdataval,
			String szmoduleval,String missing1, int missinter1, 
			boolean btakelog1, int filterduplicates) throws Exception {
		    
		//szsurvival = "D:\\Personal issues\\CREAM\\CREAM tool\\CREAM_survival v.1.0\\TCGA survival  data/TCGA-COAD.survival.tsv";
	   //szdataval = "D:\\Personal issues\\CREAM\\CREAM tool\\CREAM_survival v.1.0\\TCGA expression data/Colon Caner (COAD).txt";
	    //szmoduleval = "D:\\Personal issues\\CREAM\\CREAM tool\\CREAM_survival v.1.0\\module files (example)/COAD test file20.txt";
		
				if (szsurvival.trim().equals("") ) {
					throw new IllegalArgumentException("No survival file given!");
				} else if (!(new File(szsurvival)).exists()) {
					throw new IllegalArgumentException("The survival file '" + szsurvival+ "' cannot be found.");
				} else if (szdataval.trim().equals("")) {
					throw new IllegalArgumentException("No expression data given!");
				} else if (!(new File(szdataval)).exists()) {
					throw new IllegalArgumentException("The expression data file '" + szdataval+ "' cannot be found.");
				} else if (szmoduleval.trim().equals("")) {
					throw new IllegalArgumentException("No module file given!");
				} else if (!(new File(szmoduleval)).exists()) {
					throw new IllegalArgumentException("The module file '" + szdataval+ "' cannot be found.");
				}
				
				expressiondata = buildset(szdataval, Integer.parseInt(missing1), missinter1, btakelog1, filterduplicates);
				expressiondata.normalization_zcore();
				geneindex = new HashMap<String, Integer>();
				for(int i=0;i<expressiondata.genenames.length;i++){
					geneindex.put(expressiondata.genenames[i].toUpperCase(), i);
				}
				
				List<File> go = getAllFile(new File("Enrichment_analysis_files"));
						
				if(go != null && go.size()>0){
					pathwayData = new RegulatorBindingData(go.get(0).toString(),
							expressiondata, SZDELIM, null, null);
					pathwayData.deDuplication();
				}

				readModuleFile(szmoduleval);
				double[][] survivalfile = readsurvivalfile(szsurvival);
				
				try {
					int groupcount = (int) (expressiondata.numcols * Double.parseDouble(division));
					Survivalanalysis caseModules = new Survivalanalysis(expressiondata, pathwayData, 
							resultlist, posilist, negalist, groupcount, survivalfile, runningtext);
		        	final Survivalanalysis finalcase = caseModules;
		        	
		        	javax.swing.SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							JFrame frame1 = new JFrame("Survival Analysis Results");
							Container theDialogContainer = frame1.getContentPane();
							theDialogContainer.setBackground(Color.white);
							JTabbedPane tabbedPane = new JTabbedPane();
								CaseTable newContentPane1 = new CaseTable(frame1, 
										finalcase, finalcase.ModuleSet, expressiondata);
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
				
				/*
				String[] datasample = expressiondata.dsamplemins;
				List<Integer> sampleindex = readfile();

				String series1 = "Highly co-expressed";  
			    String series2 = "Lowly co-expressed";  
			    DefaultCategoryDataset dataset = new DefaultCategoryDataset();  
			    dataset.addValue(200, series1, "1");  
			    dataset.addValue(150, series1, "2");  
			    dataset.addValue(100, series1, "3");  
			    dataset.addValue(210, series1, "4");  
			    dataset.addValue(240, series1, "5");
			    dataset.addValue(195, series1, "6");  
			    dataset.addValue(245, series1, "7");  
			    
			    dataset.addValue(150, series2, "1");  
			    dataset.addValue(130, series2, "2");  
			    dataset.addValue(95, series2, "3");  
			    dataset.addValue(195, series2, "4");  
			    dataset.addValue(200, series2, "5");  
			    dataset.addValue(180, series2, "6");  
			    dataset.addValue(230, series2, "7"); 
				
				try {
		        	javax.swing.SwingUtilities.invokeLater(new Runnable() {
		        		public void run() {
		        			String title = "Kaplan-Meier survival plot";
		        			String pvalue = "Log-rank test p-value = 0.03";
		        			SurvivalCurve survivaltPane1 = new SurvivalCurve(title, dataset, pvalue);
							survivaltPane1.setAlwaysOnTop(true);  
							survivaltPane1.pack();  
							survivaltPane1.setSize(600, 500);
							//survivaltPane1.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);  
							survivaltPane1.setVisible(true);  
						}
					});
		            } catch (Exception e) {
						e.printStackTrace();
				}
		*/
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
            	} else {
            		System.out.println(lineTxt);
            	}
            }
            bufferedReader1.close();
    } catch (Exception e) {
    	throw new IllegalArgumentException("Read module file error.");
    }	
	}
	
	public double[][] readsurvivalfile (String filename){
		double[][] result = new double[expressiondata.numcols][2];
		List<String> sample = new ArrayList<String>();
		List<String> statue = new ArrayList<String>();
		List<String> time = new ArrayList<String>();
		try {
			BufferedReader bufferedReader1 = new BufferedReader(new FileReader(filename));   
			String lineTxt = bufferedReader1.readLine();
			while((lineTxt = bufferedReader1.readLine()) != null){
				if(Arrays.asList(lineTxt.split("\t")).size() > 2){
					sample.add(Arrays.asList(lineTxt.split("\t")).get(0));
					statue.add(Arrays.asList(lineTxt.split("\t")).get(1));
					time.add(Arrays.asList(lineTxt.split("\t")).get(2));
				}
			}
			bufferedReader1.close();
		}catch (Exception e) {
	    	throw new IllegalArgumentException("Read survival file error.");
	    }
		int matched = 0;
		for(int i=0;i<expressiondata.numcols;i++) {
			String submitter = expressiondata.dsamplemins[i];
			for(int j=0;j<sample.size();j++) {
				if(sample.get(j).equals(submitter)) {
					matched++;
					if(statue.get(j).equalsIgnoreCase("dead") || statue.get(j).equalsIgnoreCase("1")) {
						result[i][0] = 1;
						result[i][1] = Double.parseDouble(time.get(j));
					} else if(statue.get(j).equalsIgnoreCase("alive") || statue.get(j).equalsIgnoreCase("0")) {
						result[i][0] = 0;
						result[i][1] = Double.parseDouble(time.get(j));
					}
					break;
				}
			}
		}
		System.out.println(matched);
		if(matched == 0) {
			throw new IllegalArgumentException("ID mathching error between the expression data and survival file.");
		}
		return result;
	}
	/**
	 * update the input static data
	 */
	public void handlepathwaysource() {
		if (survivalField.isEditable()) {
			szuserFileField1 = survivalField.getText(); 
		}

		npathwaysourcecb = survivalsourcecb.getSelectedIndex(); 

		if (npathwaysourcecb >= 1) {
			survivalField.setText(fc
					+ System.getProperty("file.separator")
					+ staticsourceArray1[npathwaysourcecb]);
			survivalField.setEditable(false);
			survivalbutton.setEnabled(false);
		} else {
			survivalField.setText(szuserFileField1);
			survivalField.setEditable(true);
			survivalbutton.setEnabled(true);
		}
	}
	


	/**
	 * define the button methods
	 */
	public void actionPerformed(ActionEvent e) {
		Object esource = e.getSource();

		if (esource == survivalsourcecb) {
			handlepathwaysource();
		} else if (esource == survivalbutton) {
			int returnVal = fc1.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc1.getSelectedFile();
				survivalField.setText(file.getAbsolutePath());
			}
		}  else if (esource == databutton) {
			int returnVal = fc2.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc2.getSelectedFile();
				dataField.setText(file.getAbsolutePath());
			}
		} else if (esource == modulebutton) {
			int returnVal = fc3.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc3.getSelectedFile();
				moduleField.setText(file.getAbsolutePath());
			}
		}  else if (esource == Run) {
			s1 = System.currentTimeMillis();
			szsurvival = survivalField.getText();
			szdataval = dataField.getText(); 
			szmoduleval = moduleField.getText(); 
			division = jpro.getValue().toString();
			
			btakelog1 = log1.isSelected();
			missing1 = j18.getValue().toString();
			
			int temp2=comboBox1.getSelectedIndex();
			if(temp2==0){
				missinter=0;
			}else if(temp2==1){
				missinter=1;
			}else{
				missinter=2;
			}
			
			int temp3=comboBox3.getSelectedIndex();
			if(temp3==0){
				filterduplicates=0;
			}else if(temp3==1){
				filterduplicates=1;
			}
			
			

			this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			final JFrame fframe = this;

			Runnable clusterrun = new Runnable() {
				public void run() {
					Run.setEnabled(false);
					try {clusterscript
						(szsurvival, szdataval, szmoduleval, 
								missing1, missinter, btakelog1, filterduplicates);
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
