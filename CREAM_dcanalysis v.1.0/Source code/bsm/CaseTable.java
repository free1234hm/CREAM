package bsm;

import heatmapframe.GetHeatMap;
import heatmapframe.HeatMap;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSpinner;
import javax.swing.JTabbedPane;
import javax.swing.JTable;
import javax.swing.SpinnerNumberModel;

import java.awt.Color;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.event.ComponentAdapter;
import java.awt.event.ComponentEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.*;
import java.util.List;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.awt.datatransfer.*;

import javax.swing.border.TitledBorder;
import javax.swing.table.*;

import bsm.DCanalysis.Module;
import bsm.core.TableModelST;
import bsm.core.TableSorter;
import bsm.core.Util;

/**
 * Class for a table that shows enrichment of TF targets along a path
 */
public class CaseTable extends JPanel implements ActionListener {
	JFrame theframe;
	DCanalysis dcresults;
	List<Module> moduleset;
	BSM_DataSet controldata, casedata;
	JPanel mapPanel, mapPanel1, mapPanel2;
	JPanel tablePanel;
	String[] columnNames;
	String[] controlheader, caseheader;
	String[][] tabledata;
	JButton copyButton;
	JButton saveButton;
	JButton pathwayButton;
	JButton showButton;
	String showText;
	boolean showheatmap;;
	String[][] totaltable;
	HeatMap controlheatmap, caseheatmap;
	String selectedid;
	
	JLabel j2 = new JLabel("Min count: ");
	SpinnerNumberModel count = new SpinnerNumberModel(new Integer(10), new Integer(0), null, new Integer(1));
	JSpinner num1 = new JSpinner(count);
	//num1.setMinimumSize(new Dimension(20,20));
	JLabel j3 = new JLabel("Min percentage:");
	SpinnerNumberModel per = new SpinnerNumberModel(new Integer(10), new Integer(0), null, new Integer(1));
	JSpinner num2 = new JSpinner(per);
	//num2.setMinimumSize(new Dimension(20,20));
	JLabel j4 = new JLabel("%");
	JScrollPane scrollPane;
	JScrollPane scrollcontrol, scrollcase;
	TableSorter sorter;
	TableSorter sortercontrol, sortercase;
	JTable table;
	JTable tablecontrol, tablecase;
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	NumberFormat nf;
	NumberFormat nf2;
	boolean bsplit;
	DecimalFormat df1;
	DecimalFormat df2;
	
	/**
	 * Constructor - builds the table
	 */
	public CaseTable(JFrame frame1, DCanalysis dcresults, 
			List<Module> moduleset, BSM_DataSet controldata, BSM_DataSet casedata) {
		
		this.theframe = frame1;
		this.dcresults = dcresults;
		this.moduleset = moduleset;
		this.controldata = controldata;
		this.casedata = casedata;
		df1 = new DecimalFormat("#0.00");
		df2 = new DecimalFormat("#0.000");
		
		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);

		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);
		
		tabledata = new String[moduleset.size()][7];	
		columnNames = new String[7];
		columnNames[0] = "Module ID";
		columnNames[1] = "Total genes";
		columnNames[2] = "Positive genes";
		columnNames[3] = "Negative genes";
		columnNames[4] = "Control Avg. correlation";
		columnNames[5] = "Case Avg. correlation";
		columnNames[6] = "DC p-value";
		for(int i=0;i<moduleset.size();i++) {
			Module module = moduleset.get(i);
			tabledata[i][0] = module.id;
			tabledata[i][1] = module.list.size()+"";
			tabledata[i][2] = module.posilist.size()+"";
			tabledata[i][3] = module.negalist.size()+"";
			tabledata[i][4] = df2.format(module.dcpvalue[0]);
			tabledata[i][5] = df2.format(module.dcpvalue[1]);
			tabledata[i][6] = module.dcpvalue[2]+"";
		}
		
		sorter = new TableSorter(new TableModelST(tabledata, columnNames));
		table = new JTable(sorter);
		sorter.setTableHeader(table.getTableHeader());
		//table.setPreferredScrollableViewportSize(new Dimension(0, 
		//Math.min((table.getRowHeight() + table.getRowMargin())* table.getRowCount(), 400)));

		//table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
		TableColumn column = table.getColumnModel().getColumn(0);
		column.setMinWidth(150);
		for(int i=1;i<columnNames.length-1;i++) {
			column = table.getColumnModel().getColumn(i);
			column.setMinWidth(50);
		}
		column = table.getColumnModel().getColumn(columnNames.length-1);
		column.setMinWidth(150);
		
		scrollPane = new JScrollPane(table);
		scrollPane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
		scrollPane.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED);
		
		tablePanel = new JPanel();
		tablePanel.setBorder(new TitledBorder(null,"Module Information",TitledBorder.LEFT,TitledBorder.TOP));
		BoxLayout layout = new BoxLayout(tablePanel, BoxLayout.Y_AXIS);
		tablePanel.setLayout(layout);
		tablePanel.add(scrollPane);
		
		addheatmap();
		add(tablePanel);
		addBottom();
		
		this.addComponentListener(new ComponentAdapter(){  //dynamic column size;
			public void componentResized(ComponentEvent e){
				for(int i=0;i<columnNames.length;i++) {
					TableColumn column = table.getColumnModel().getColumn(i);
					column.setPreferredWidth(Math.max(50, theframe.getWidth()/columnNames.length));
				}
		  }
	   });
		
		table.addMouseListener(new MouseAdapter() {
		     public void mouseClicked(MouseEvent e) {
		    	 if(table.getSelectedRow() != -1){
		    		 mapPanel1.removeAll();
		    		 mapPanel1.repaint();
		    		 mapPanel2.removeAll();
		    		 mapPanel2.repaint();
		    		 
		    		 selectedid = (String) table.getValueAt(table.getSelectedRow(), 0);
		    		 Module selectmodule = null;
		    		 for(int i=0;i<moduleset.size();i++) {
		    			 if(selectedid.equals(moduleset.get(i).id)) {
		    				 selectmodule = moduleset.get(i);
		    			 }
		    		 }
		    		 
		    		 controlheader = new String[controldata.numcols+1];
		    		 controlheader[0] = "ID";
	    			 for(int i=0;i<controldata.numcols;i++) controlheader[i+1] = controldata.dsamplemins[i];
	    			 caseheader = new String[casedata.numcols+1];
	    			 caseheader[0] = "ID";
	    			 for(int i=0;i<casedata.numcols;i++) caseheader[i+1] = casedata.dsamplemins[i];
		    		 
		    		 GetHeatMap heatmap = new GetHeatMap();
		    		 
		    		 if(selectmodule != null && selectmodule.controltable != null && selectmodule.casetable != null) {	 
		    			 double[][] normcontrol = new double[selectmodule.controltable.length][controldata.numcols];
		    			 for(int i=0;i<selectmodule.controltable.length;i++) {
		    				 double largest = Double.MIN_VALUE;
						     double smallest = Double.MAX_VALUE;
		    				 for(int j=1;j<selectmodule.controltable[0].length;j++){
		    					 largest = Math.max(Double.parseDouble(selectmodule.controltable[i][j]), largest);
								 smallest = Math.min(Double.parseDouble(selectmodule.controltable[i][j]), smallest);
		    				 }
		    				 double range = largest - smallest;
		    				 for(int j=0;j<controldata.numcols;j++) normcontrol[i][j] = 
		    						 (Double.parseDouble(selectmodule.controltable[i][j+1])-smallest)/range;
				    	}
		    			 double[][] normcase = new double[selectmodule.casetable.length][casedata.numcols];
		    			 for(int i=0;i<selectmodule.casetable.length;i++) {
		    				 double largest = Double.MIN_VALUE;
						     double smallest = Double.MAX_VALUE;
		    				 for(int j=1;j<selectmodule.casetable[0].length;j++){
		    					 largest = Math.max(Double.parseDouble(selectmodule.casetable[i][j]), largest);
								 smallest = Math.min(Double.parseDouble(selectmodule.casetable[i][j]), smallest);
		    				 }
		    				 double range = largest - smallest;
		    				 for(int j=0;j<casedata.numcols;j++) normcase[i][j] = 
		    						 (Double.parseDouble(selectmodule.casetable[i][j+1])-smallest)/range;
		    				 } 
		    			 try {
		    				 controlheatmap = heatmap.heatmap(normcontrol);
		    				 caseheatmap = heatmap.heatmap(normcase);
		    				 } catch (Exception e2) {
		    					 e2.printStackTrace();
		    				}
		    			 
		    			    sortercontrol = new TableSorter(new TableModelST(selectmodule.controltable, controlheader));
			    			tablecontrol = new JTable(sortercontrol);
			    			sortercontrol.setTableHeader(tablecontrol.getTableHeader());
			    			tablecontrol.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			    			TableColumn column;
			    			for(int i=0;i<controlheader.length;i++){
			    				column = tablecontrol.getColumnModel().getColumn(i);
			    				column.setMaxWidth(120);
			    			}
			    			scrollcontrol = new JScrollPane(tablecontrol);
			    			scrollcontrol.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
			    			scrollcontrol.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
			    			
			    			sortercase = new TableSorter(new TableModelST(selectmodule.casetable, caseheader));
			    			tablecase = new JTable(sortercase);
			    			sortercase.setTableHeader(tablecase.getTableHeader());
			    			tablecase.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
			    			for(int i=0;i<caseheader.length;i++){
			    				column = tablecase.getColumnModel().getColumn(i);
			    				column.setMaxWidth(120);
			    			}
			    			scrollcase = new JScrollPane(tablecase);
			    			scrollcase.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
			    			scrollcase.setVerticalScrollBarPolicy(JScrollPane.VERTICAL_SCROLLBAR_ALWAYS);
			    			
			       			if(showheatmap){
								mapPanel1.add(controlheatmap);
								mapPanel2.add(caseheatmap);
				    		 }else{
				    			mapPanel1.add(scrollcontrol);
				    			mapPanel2.add(scrollcase);
				    		 }
		    		 }
		    		 mapPanel1.updateUI();
		    		 mapPanel2.updateUI();
		    		 //mapPanel.updateUI();
		    	 }
		     }
		});
	}

  private void addheatmap() {
	  
	  mapPanel = new JPanel();
	  mapPanel.setBorder(new TitledBorder(null,"Gene Table",TitledBorder.LEFT,TitledBorder.TOP));
	  BoxLayout layout = new BoxLayout(mapPanel, BoxLayout.X_AXIS);
	  mapPanel.setLayout(layout);
	  mapPanel.setPreferredSize(new Dimension(0, 400));
	  add(mapPanel);
	  
      mapPanel1 = new JPanel();
	  mapPanel1.setBorder(new TitledBorder(null,"Control expression matrix",TitledBorder.CENTER,TitledBorder.TOP));
	  BoxLayout layout1 = new BoxLayout(mapPanel1, BoxLayout.X_AXIS);
	  mapPanel1.setLayout(layout1);
	  mapPanel1.setPreferredSize(new Dimension(0, 400));
	  mapPanel.add(mapPanel1);
	  
	  mapPanel2 = new JPanel();
	  mapPanel2.setBorder(new TitledBorder(null,"Case expression matrix",TitledBorder.CENTER,TitledBorder.TOP));
	  BoxLayout layout2 = new BoxLayout(mapPanel2, BoxLayout.X_AXIS);
	  mapPanel2.setLayout(layout2);
	  mapPanel2.setPreferredSize(new Dimension(0, 400));
	  mapPanel.add(mapPanel2);

	}
	
  private double[][] transpose(double[][] mm){
	  double[][] result = new double[mm[0].length][mm.length];
	  for(int i=0;i<mm.length;i++) {
		  for(int j=0;j<mm[0].length;j++) {
			  result[j][i] = mm[i][j];
		  }
	  }
	  return result;
  }
  
	/**
	 * Helper function that adds information displayed at the bottom of the
	 * table information window
	 */
  
	private void addBottom() {

		copyButton = new JButton("Copy Table", Util.createImageIcon("Copy16.gif"));
		copyButton.setActionCommand("copy");
		copyButton.setMinimumSize(new Dimension(800, 20));
		copyButton.addActionListener(this);

		saveButton = new JButton("Save Table", Util.createImageIcon("Save16.gif"));
		saveButton.setActionCommand("save");
		saveButton.setMinimumSize(new Dimension(800, 20));
		saveButton.addActionListener(this);
		
		pathwayButton = new JButton("Enrichment Analysis");
		pathwayButton.setActionCommand("enrichment");
		pathwayButton.setMinimumSize(new Dimension(800, 20));
		pathwayButton.addActionListener(this);
		
		showheatmap = false;
		if (showheatmap) {
			showText = "HeatMap";
		} else {
			showText = "GeneList";
		}
		showButton = new JButton(showText);
		showButton.setActionCommand("show");
		showButton.setMinimumSize(new Dimension(800, 20));
		showButton.addActionListener(this);
		
		JPanel buttonPanel = new JPanel();
		//buttonPanel.setBackground(Color.white);
		buttonPanel.add(showButton);
		buttonPanel.add(pathwayButton);
		buttonPanel.add(saveButton);
		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE, 20));
		add(buttonPanel);
	}

	/**
	 * Writes the content of the table to a file specified through pw
	 */
	public void printFile(PrintWriter pw) {

		for (int ncol = 0; ncol < columnNames.length - 1; ncol++) {
			pw.print(columnNames[ncol] + "\t");
		}
		pw.println(columnNames[columnNames.length - 1]);

		for (int nrow = 0; nrow < tabledata.length; nrow++) {
			for (int ncol = 0; ncol < tabledata[nrow].length - 1; ncol++) {
				pw.print(sorter.getValueAt(nrow, ncol) + "\t");
			}
			pw.println(sorter.getValueAt(nrow, columnNames.length - 1));
		}
	}
	
	/**
	 * Copies the content of the table to the clipboard
	 */
	public void writeToClipboard() {
		StringBuffer sbuf = new StringBuffer();
		for (int ncol = 0; ncol < columnNames.length - 1; ncol++) {
			sbuf.append(columnNames[ncol] + "\t");
		}
		sbuf.append(columnNames[columnNames.length - 1] + "\n");

		for (int nrow = 0; nrow < tabledata.length; nrow++) {
			for (int ncol = 0; ncol < tabledata[nrow].length - 1; ncol++) {
				sbuf.append(sorter.getValueAt(nrow, ncol) + "\t");
			}
			sbuf.append(sorter.getValueAt(nrow, columnNames.length - 1) + "\n");
		}
		// get the system clipboard
		Clipboard systemClipboard = Toolkit.getDefaultToolkit()
				.getSystemClipboard();
		// set the textual content on the clipboard to our
		// Transferable object

		Transferable transferableText = new StringSelection(sbuf.toString());
		systemClipboard.setContents(transferableText, null);
	}
	
	/**
	 * Responds to buttons being pressed on the interface
	 */
	public void actionPerformed(ActionEvent e) {
		
		String szCommand = e.getActionCommand();
		if (szCommand.equals("copy")) {
			writeToClipboard();
		} else if (szCommand.equals("show")) {
			
			if (showheatmap) {
				showheatmap = false;
				if (showheatmap) {
					showText = "HeatMap";
				} else {
					showText = "GeneList";
				}
				showButton.setText(showText);
			} else {
				showheatmap = true;
				if (showheatmap) {
					showText = "HeatMap";
				} else {
					showText = "GeneList";
				}
				showButton.setText(showText);
			}
			Showheatmap();
			
		} else if (szCommand.equals("save")) {
			try {
				int nreturnVal = Main_interface.theChooser.showSaveDialog(this);
				if (nreturnVal == JFileChooser.APPROVE_OPTION) {
					File f = Main_interface.theChooser.getSelectedFile();
					PrintWriter pw = new PrintWriter(new FileOutputStream(f));
					if (szCommand.equals("save")) {
						printFile(pw);
					}
					pw.close();
				}
			} catch (final FileNotFoundException fex) {
				javax.swing.SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JOptionPane.showMessageDialog(null, fex.getMessage(),
								"Exception thrown", JOptionPane.ERROR_MESSAGE);
					}
				});
				fex.printStackTrace(System.out);
			}
		} else if (szCommand.equals("enrichment")) {
			selectedid = (String) table.getValueAt(table.getSelectedRow(), 0);
			Module selectmodule = null;
			for(int i=0;i<moduleset.size();i++) {
				if(selectedid.equals(moduleset.get(i).id)) {
					selectmodule = moduleset.get(i);
					}
				}

				String[][] pathwaylist = selectmodule.pathwayenrichment;
				javax.swing.SwingUtilities.invokeLater(new Runnable() {
					public void run() {
						JFrame pathwayframe = new JFrame("Enrichment analysis");
						JDialog frame = new JDialog(pathwayframe, "Enrichment analysis", true);
						frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
						frame.setLocation(20, 50);
						Container theDialogContainer = frame.getContentPane();
						theDialogContainer.setBackground(Color.white);
						JTabbedPane tabbedPane = new JTabbedPane();
						
						if(pathwaylist!=null && pathwaylist.length>0) {
							TSMinerGui_PathwayTable newContentPane1 = new TSMinerGui_PathwayTable(
									pathwaylist, frame, pathwayframe);
							newContentPane1.setOpaque(true); // content panes must be opaque
							tabbedPane.addTab("pathway", null, newContentPane1,"pathway");
							theDialogContainer.add(tabbedPane);
						}

						pathwayframe.setContentPane(theDialogContainer);
						pathwayframe.pack();
						pathwayframe.setVisible(true);
					}
				});
			
		}
	}

	/**
	 * Converts the value of dval to a String that is displayed on the table
	 */
	public static String doubleToSz(double dval) {
		String szexp;
		double dtempval = dval;
		int nexp = 0;

		NumberFormat nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(3);
		nf2.setMaximumFractionDigits(3);

		NumberFormat nf1 = NumberFormat.getInstance(Locale.ENGLISH);
		nf1.setMinimumFractionDigits(2);
		nf1.setMaximumFractionDigits(2);

		if (dval <= 0) {
			szexp = "0.000";
		} else {
			while ((dtempval < 0.9995) && (dtempval > 0)) {
				nexp--;
				dtempval = dtempval * 10;
			}

			if (nexp < -2) {
				dtempval = Math.pow(10, Math.log(dval) / Math.log(10) - nexp);
				szexp = nf1.format(dtempval) + "e" + nexp;

			} else {
				szexp = nf2.format(dval);
			}
		}

		return szexp;
	}
	
	class Pathway_pval {
		String pathwayID;
		double pvalue;
		Pathway_pval(String pathwayID, double pvalue) {
			this.pathwayID = pathwayID;
			this.pvalue = pvalue;
		}
	}


    
    public String[][] BubbleSort_increase(String[][] r, Integer n, Integer col) //����ð������
 	{
 		 int low = 0;   
 		    int high= n - 1; //���ñ����ĳ�ʼֵ  
 		    String[] tmp;
 		    int j;  
 		    while (low < high) {  
 		        for (j= low; j< high; ++j) //����ð��,�ҵ������  
 		            if (Double.parseDouble(r[j][col])> Double.parseDouble(r[j+1][col])) {  
 		                tmp = r[j]; r[j]=r[j+1];r[j+1]=tmp;  
 		            }   
 		        --high;                 //�޸�highֵ, ǰ��һλ  
 		        for ( j=high; j>low; --j) //����ð��,�ҵ���С��  
 		            if (Double.parseDouble(r[j][col])<Double.parseDouble(r[j-1][col])) {  
 		                tmp = r[j]; r[j]=r[j-1];r[j-1]=tmp;  
 		            }  
 		        ++low; //�޸�lowֵ,����һλ  
 		    }   
 	return r;
 	}
    
    public String[][] BubbleSort_decrease(String[][] r, Integer n, Integer col) //����ð������
 	{
 		 int low = 0;   
 		    int high= n - 1; //���ñ����ĳ�ʼֵ  
 		    String[] tmp;
 		    int j;  
 		    while (low < high) {  
 		        for (j= low; j< high; ++j) //����ð��,�ҵ������  
 		            if (Double.parseDouble(r[j][col])< Double.parseDouble(r[j+1][col])) {  
 		                tmp = r[j]; r[j]=r[j+1];r[j+1]=tmp;  
 		            }   
 		        --high;                 //�޸�highֵ, ǰ��һλ  
 		        for ( j=high; j>low; --j) //����ð��,�ҵ���С��  
 		            if (Double.parseDouble(r[j][col])>Double.parseDouble(r[j-1][col])) {  
 		                tmp = r[j]; r[j]=r[j-1];r[j-1]=tmp;  
 		            }  
 		        ++low; //�޸�lowֵ,����һλ  
 		    }   
 	return r;
 	}

    public Integer getrank(List<Double> pvalue, double value){
    	int cc = 1;
    	for(int i=0;i<pvalue.size();i++){
    		if(pvalue.get(i)<value){
    			cc++;
    		}
    	}
    	return cc;
    }
    private void Showheatmap() {
    	if(table.getSelectedRow() != -1){
        	 mapPanel1.removeAll();
       		 mapPanel1.repaint();
       		 mapPanel2.removeAll();
      		 mapPanel2.repaint();
       		 if(showheatmap){
    			mapPanel1.add(controlheatmap);
				mapPanel2.add(caseheatmap);
       		 }else{
       			mapPanel1.add(scrollcontrol);
       			mapPanel1.add(scrollcase);
       		 }
       		mapPanel1.updateUI();
       		mapPanel2.updateUI();
       		//mapPanel.updateUI();
         }
	}
    
}