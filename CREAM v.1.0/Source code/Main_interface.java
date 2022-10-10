package bsm;


import bsm.core.*;

import java.awt.*;
import java.awt.event.*;
import java.io.*;

import javax.swing.*;
import javax.swing.border.TitledBorder;

import java.util.*;
import java.util.List;

/**
 * Class implementing the main input interface
 */
public class Main_interface extends JFrame implements ActionListener {
	
	static final String SZDELIM = "|;,";
	static final boolean BDEBUG = false;
	
	BSM_DataSet theDataSet;
	List<RegulatorBindingData> KnownModules;
	int numbits;
	
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
	ListDialog theRepeatList;
	static Color lightBlue = new Color(190, 255, 190);
	Vector<String> knownmodules;

	// Search Options
	static double dCONVERGENCEDEF = 0.01;
	static double dMinScoreDEF = 0.0;
	static double dDELAYPATHDEF = .15;
	static double dDMERGEPATHDEF = .15;
	static double dPRUNEPATHDEF = .15;
	static double dTHRESHOLD = 0.000001;
	static double dPvalueshold = 0.01;

	// Filtering
	JRadioButton log1, log2, go1, go2;
	static int nnormalizeDEF = 1;
	static int ngoDEF = 0;
	static int ntfDEF = 1;
	static int nmirnaDEF = 1;
	ButtonGroup normGroup = new ButtonGroup();
	ButtonGroup goGroup = new ButtonGroup();
	ButtonGroup tfGroup = new ButtonGroup();
	ButtonGroup mirnaGroup = new ButtonGroup();
	
	static boolean busego;
	
	String labels1[] = {" Min value", " Mean value", " Zero"};
	JComboBox comboBox1 = new JComboBox(labels1);
	String labels3[] = {" Max value", " Mean value"};
	JComboBox comboBox3 = new JComboBox(labels3);
	//String labels4[] = {" Yes"," No"};
	//JComboBox comboBox4 = new JComboBox(labels4);
	JLabel filterlable;
	static int missinter = 0; // Set the missing value as
	static int filterduplicates = 0; // Set the duplicate genes as
	static String szPrefilteredDEF = "";
	static int nMaxMissingDEF = 0;
	static double dMinExpressionDEF = 1;
	static double dMinCorrelationRepeatsDEF = 0;

	// GUI
	static long s1;
	static boolean bendsearch = false;
	static String szexpressionval;
	static String szsavedval;
	static List<List<String>> savedmodel;
	static String missing;
	static boolean btakelog;
	static String epsilon;
	static String szmaxchild;
	static String szminchild;
	static String mincorr;
	static String DEthreshold;
	static String FDRthreshold;
	static String FCthreshold;

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
	static Color defaultColor;
	static String[] staticsourceArray1 = { "User provided" };
	static String[] staticsourceArray2 = { "User provided" };
	static String[] staticsourceArray3 = { "User provided" };

	JTextField goField;
	JTextField extraField;
	JTextField xrefField;
	JTextField categoryIDField;
	JTextField taxonField;
	JTextField evidenceField;
	JButton infoButton = new JButton(Util.createImageIcon("About16.gif"));
	static JFileChooser fc1 = new JFileChooser(new File("Known gene modules"));
	static JFileChooser fc2 = new JFileChooser(new File("Expression data (sample)"));
	static JFileChooser fc3 = new JFileChooser(new File("Saved module file (sample)"));
	

	static Container contentPane;
	JPanel textpanel;
	JPanel panel4;

	JButton pathwaybutton= new JButton();
	JButton timeseriesbutton= new JButton();
	JButton savedmodelbutton= new JButton();
	JTextArea runningtext;
	JButton Run = new JButton();
	JButton currentButton = new JButton();
	JButton endSearchButton = new JButton();
	JTextField pathwayField;
	JTextField expressionField;
	JTextField savedmodelField;
	final JSpinner j18, j19, j20, j21, j22;

	
	public Main_interface() throws FileNotFoundException, IOException {
		super("CREAM v.1.0");

		contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.Y_AXIS);
		contentPane.setLayout(layout);

		JPanel panel1=new JPanel();
		panel1.setBorder(new TitledBorder(null,"Load Data",TitledBorder.LEFT,TitledBorder.TOP));
		panel1.setLayout(new BorderLayout());
		
/*****************************************************************/
		
		JLabel j3=new JLabel("Load known modules:");
		j3.setBounds(15,20,300,35);
		pathwayField=new JTextField(szoptionalFile, JLabel.TRAILING);
		pathwayField.setBounds(170,23,400,25);
		pathwayField.setBorder(BorderFactory.createLineBorder(Color.black));
		pathwayField.setOpaque(true);
		pathwayField.setBackground(Color.white);
		pathwaybutton.setText("Load File");
		pathwaybutton.setHideActionText(true);
		pathwaybutton.addActionListener(this);
		pathwaybutton.setBounds(575,23,90,25);
		
		/*****************************************************************/
		
		JLabel j7=new JLabel("Load expression file:");
		j7.setBounds(15,55,300,35);
		expressionField=new JTextField(szDataFileDEF, JLabel.TRAILING);
		expressionField.setBounds(170,58,400,25);
		expressionField.setBorder(BorderFactory.createLineBorder(Color.black));
		expressionField.setOpaque(true);
		expressionField.setBackground(Color.white);
		timeseriesbutton.setText("Load File");
		timeseriesbutton.setHideActionText(true);
		timeseriesbutton.addActionListener(this);
		timeseriesbutton.setBounds(575,58,90,25);
		
		/*****************************************************************/
		
		JLabel j8=new JLabel("Load saved model file:");
		j8.setBounds(15,90,300,35);
		savedmodelField=new JTextField(szoptionalFile, JLabel.TRAILING);
		savedmodelField.setBounds(170,93,400,25);
		savedmodelField.setBorder(BorderFactory.createLineBorder(Color.black));
		savedmodelField.setOpaque(true);
		savedmodelField.setBackground(Color.white);
		savedmodelbutton.setText("Load File");
		savedmodelbutton.setHideActionText(true);
		savedmodelbutton.addActionListener(this);
		savedmodelbutton.setBounds(575,93,90,25);
		
		/*****************************************************************/

		panel1.setLayout(null);
		//panel1.add(j11); 
		//panel1.add(go1);
		//panel1.add(go2);
		goGroup.add(go1);
		goGroup.add(go2);

		
		panel1.add(j3);
		panel1.add(j7);
		panel1.add(j8);
		
		panel1.add(pathwayField);
		panel1.add(expressionField);
		panel1.add(savedmodelField);
		
		panel1.add(pathwaybutton);
		panel1.add(timeseriesbutton);
		panel1.add(savedmodelbutton);
		
        JScrollPane   scrollpanel1   =   new   JScrollPane(panel1, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel1.setBounds(100, 100, 745, 70);
		panel1.setPreferredSize(new Dimension(scrollpanel1.getWidth() - 50, scrollpanel1.getHeight()*2));
		
		/*******************************************************/
		
		JPanel panel3=new JPanel();
		JScrollPane   scrollpanel2   =   new   JScrollPane(panel3, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollpanel2.setBounds(100, 100, 745, 100);
        panel3.setPreferredSize(new Dimension(scrollpanel2.getWidth() - 50, scrollpanel2.getHeight()*2));
		panel3.setBorder(new TitledBorder(null,"Set Parameters",TitledBorder.LEFT,TitledBorder.TOP));
		panel3.setLayout(new BorderLayout());
		
		JPanel panel31=new JPanel();
		panel31.setBorder(new TitledBorder(null,"Data Preprocessing",TitledBorder.LEFT,TitledBorder.TOP));
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
		panel32.setBorder(new TitledBorder(null,"Module Detecting",TitledBorder.LEFT,TitledBorder.TOP));
		panel32.setBounds(370,20,365,175);
		panel32.setLayout(null);

		JLabel minir=new JLabel("Correlation threshold:");
		minir.setBounds(15,20,300,35);
		SpinnerNumberModel minipearson = new SpinnerNumberModel(0.4, -1, 1, 0.1);
		j22 = new JSpinner(minipearson);
		j22.setBounds(150,25,50,22);
		
		JLabel knn=new JLabel("Min improvement of likelihood:");
		knn.setBounds(15,55,200,35);
		SpinnerNumberModel con = new SpinnerNumberModel(0.005, 0, 0.01, 0.001);
		j20 = new JSpinner(con);
		j20.setBounds(200,60,50,22);
		
		JLabel nummodule1=new JLabel("Max");
		nummodule1.setBounds(15,90,30,35);
		SpinnerNumberModel converenge = new SpinnerNumberModel(numchildDEF, 2, null, 1);
		j21 = new JSpinner(converenge);
		j21.setBounds(42,95,50,22);
		JLabel nummodule2=new JLabel("and min");
		nummodule2.setBounds(95,90,50,35);	
		SpinnerNumberModel minigene = new SpinnerNumberModel(1, 1, null, 1);
		j19 = new JSpinner(minigene);
		j19.setBounds(145,95,50,22);	
		JLabel nummodule3=new JLabel("number of modules");
		nummodule3.setBounds(198,90,200,35);
		
		JLabel j11=new JLabel("Use known modules for gene assignment:");
		j11.setBounds(15,125,300,35);
		go1 = new JRadioButton("Yes");
		go2 = new JRadioButton("No");
		go1.setBounds(260, 131, 50, 22);
		go2.setBounds(310, 131, 50, 22);
		if (ngoDEF == 0) {
			go1.setSelected(true);
		} else if (ngoDEF == 1) {
			go2.setSelected(true);
		}

		panel32.add(knn);
		panel32.add(nummodule1);
		panel32.add(nummodule2);
		panel32.add(nummodule3);
		panel32.add(minir);
		panel32.add(j19);
		panel32.add(j20);
		panel32.add(j21);
		panel32.add(j22);
		panel32.add(j11);
		panel32.add(go1);
		panel32.add(go2);
		goGroup.add(go1);
		goGroup.add(go2);

		panel3.setLayout(null);
		panel3.add(panel31); 
		panel3.add(panel32);
		/********************************************************/
		textpanel=new JPanel();
		textpanel.setPreferredSize(new Dimension(745, 150));
		textpanel.setBorder(new TitledBorder(null,"Search Gene Modules",TitledBorder.LEFT,TitledBorder.TOP));
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
		endSearchButton.setText("Stop Searching");
		endSearchButton.setEnabled(false);
		endSearchButton.addActionListener(this);
		infoButton.addActionListener(this);
		infoButton.setBackground(gray);
		
		JPanel panel6=new JPanel();
		panel6.setPreferredSize(new Dimension(745, 50));
		panel6.add(Run);
		panel6.add(endSearchButton);
		panel6.add(infoButton);
		contentPane.add(scrollpanel1);
		contentPane.add(scrollpanel2);
		contentPane.add(textpanel);
		contentPane.add(panel6);
		
		theRepeatList = new ListDialog(this, Main_interface.vRepeatFilesDEF,
				pathwaybutton, Main_interface.fc1);
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
			Integer maxcase, Integer missinter, boolean btakelog) 
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
	public void clusterscript(Vector<String> knownmodules, String szexp1val, String szsavedval,
			String missing, int missinter, boolean btakelog) throws Exception {
		
		//String aa = "D:/Personal issues/CREAM/data_annotation/GO&pathway/KEGG pathways/KEGG_pathway_geneset_hsa.txt";
		//knownmodules.add(aa);
		//szexp1val = "D:\\Personal issues\\CREAM\\data_expression\\TCGA\\GDC\\colon cancer\\pathologic_TNM/COAD_T4.txt";
		//szsavedval = "D:\\Personal issues\\CREAM\\result\\GDC\\epsilon=0.005\\1_Module files\\CREAM\\COAD\\T4/finalmodel.txt";
		
				if (szexp1val.trim().equals("")) {
					throw new IllegalArgumentException("No expression data file given!");
				} else if (!(new File(szexp1val)).exists()) {
					throw new IllegalArgumentException("The expression data file '" + szexp1val+ "' cannot be found.");
				}
				
				bendsearch = false;			
				theDataSet = buildset(szexp1val, Integer.parseInt(missing), missinter, btakelog);
				theDataSet.normalization();
				
				if(knownmodules.size()>0) {
					KnownModules = new ArrayList<RegulatorBindingData>();
					try {
						for(String file:knownmodules) {
							RegulatorBindingData pathwayData = new RegulatorBindingData(file,
									theDataSet, SZDELIM, null, null);
							pathwayData.deDuplication();
							KnownModules.add(pathwayData);
						}
					}catch (Exception e) {
						throw new IllegalArgumentException("Read known module file error, please check the knwon module files.");
				    }
				}
				
				if(szsavedval != null && szsavedval.length()>0 && !szsavedval.equals(" Optional")){
					try {
						savedmodel = readFile(szsavedval);
					} catch (Exception e) {
						throw new IllegalArgumentException("Please check the saved case model file.");
				    }
				}
				
				try {
		        	//System.out.println(savedcase.size()+" "+savedcase.get(0).size());
		        	LearningCase caseModules = new LearningCase(theDataSet, KnownModules, 
		        			busego, savedmodel, szmaxchild, szminchild, mincorr,
		        			epsilon, runningtext, endSearchButton);
		        	
		        	final LearningCase finalcase = caseModules;
		        	
		        	javax.swing.SwingUtilities.invokeLater(new Runnable() {
						public void run() {
							JFrame frame1 = new JFrame("Resulting Modules");
							Container theDialogContainer = frame1.getContentPane();
							theDialogContainer.setBackground(Color.white);
							JTabbedPane tabbedPane = new JTabbedPane();
								CaseTable newContentPane1 = new CaseTable(frame1,
										finalcase, finalcase.treeptr);
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

	/**
	 * define the button methods
	 */
	public void actionPerformed(ActionEvent e) {
		Object esource = e.getSource();

		if (esource == endSearchButton) {
			bendsearch = true;
			runningtext.append("End Search Requested. Search Will End Soon..."+"\n");
			runningtext.paintImmediately(runningtext.getBounds());
			endSearchButton.setEnabled(false);
		} else if (esource == pathwaybutton) {
			theRepeatList.setLocation(this.getX() + 75, this.getY() + 100);
			theRepeatList.setVisible(true);
		} else if (esource == timeseriesbutton) {
			int returnVal = fc2.showOpenDialog(this);
			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc2.getSelectedFile();
				expressionField.setText(file.getAbsolutePath());
			}
		} else if (esource == savedmodelbutton) {
			int returnVal = fc3.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc3.getSelectedFile();
				savedmodelField.setText(file.getAbsolutePath());
			}
		}  else if (esource == Run) {
			s1 = System.currentTimeMillis();
			szexpressionval = expressionField.getText(); 
			knownmodules = theRepeatList.data;
			szsavedval = savedmodelField.getText(); 
			busego = go1.isSelected();
			missing = j18.getValue().toString();
			epsilon = j20.getValue().toString();
			szmaxchild = j21.getValue().toString();
			szminchild = j19.getValue().toString();
			mincorr = j22.getValue().toString();
			
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
			
			btakelog = log1.isSelected();

			this.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));

			final JFrame fframe = this;

			Runnable clusterrun = new Runnable() {
				public void run() {
					Run.setEnabled(false);
					try {
						clusterscript(knownmodules, szexpressionval, 
								szsavedval, missing, missinter, btakelog);
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
