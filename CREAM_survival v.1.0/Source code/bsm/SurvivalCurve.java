package bsm;


import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.data.category.DefaultCategoryDataset;


/**
 * Class for a table that shows enrichment of TF targets along a path
 */
public class SurvivalCurve extends JFrame {

	/**
	 * Constructor - builds the table
	 */
	public SurvivalCurve(String title, DefaultCategoryDataset dataset, String pvalue) {
		super(title);
	    JFreeChart chart = ChartFactory.createLineChart(  
	            pvalue, // Chart title  
	            "Time", // X-Axis Label  
	            "Survival probability (%)", // Y-Axis Label  
	            dataset  
	            );
	    ChartPanel panel = new ChartPanel(chart);
	    setContentPane(panel); 
	}
}