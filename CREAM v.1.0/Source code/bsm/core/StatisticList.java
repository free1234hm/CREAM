package bsm.core;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

import javax.swing.filechooser.*;

import java.awt.image.*;

public class StatisticList extends JDialog implements ActionListener {


	String labels[] = { "Student's t-test", "Wilcoxon rank sum test"};
	JComboBox comboBox3 = new JComboBox(labels);
	JButton okButton;
	Color defaultColor;
	public static Color buttonColor = new Color(255, 246, 143);
	JFrame frame;
//	Vector<String> names = new Vector<String>();
	Color bColor;

	/**
	 * Constructors the Repeat interface
	 */
	public StatisticList(JFrame frame, Color bColor, Color defaultColor) {
		super(frame, "Statistic Methods", false);
		this.frame = frame;
		this.bColor = bColor;
		this.defaultColor = defaultColor;

		
		JLabel label = new JLabel("Statistic Method:");
		comboBox3.setPreferredSize(new Dimension(200, 28));
		comboBox3.setMaximumSize(new Dimension(200, 28));
		okButton = new JButton("OK");
		okButton.setActionCommand("ok");
		okButton.addActionListener(this);
		okButton.setPreferredSize(new Dimension(100, 28));
		okButton.setMaximumSize(new Dimension(100, 28));
		
		JPanel listPane1 = new JPanel();
		listPane1.add(label);
		listPane1.add(comboBox3);
		listPane1.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		
		JPanel listPane2 = new JPanel();
		listPane2.add(okButton);
		listPane2.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

		Container contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.LINE_AXIS);
		contentPane.setLayout(layout);
		contentPane.setBackground(bColor);
		contentPane.add(listPane1);
		contentPane.add(listPane2);
		pack();
	}

	/**
	 * Handle clicks on the Set and Cancel buttons.
	 */
	public void actionPerformed(ActionEvent e) {
		String szcommand = e.getActionCommand();
		if (szcommand.equals("ok")) {
			setVisible(false);
		}
	}

}
