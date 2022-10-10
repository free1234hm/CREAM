package bsm.core;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.io.*;

import javax.swing.filechooser.*;

import java.awt.image.*;

/*
 * Classed used for the repeat dialog interface. It is modeled after this sample code:
 * http://java.sun.com/docs/books/tutorial/uiswing/examples/components/ListDialogRunnerProject/src/components/ListDialog.java
 * here is the copyright for that code
 *
 * Copyright (c) 1995 - 2008 Sun Microsystems, Inc.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *
 *   - Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *
 *   - Neither the name of Sun Microsystems nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
public class ListDialog extends JDialog implements ActionListener {

//	private static ListDialog dialog;
//	private static String value = "";
	private JList list;
	private JList list2;
	public Vector<String> data;
	public Vector<String> statisticlist;
	JButton mainButton;
	JButton removeButton;
	JButton viewButton;
	JButton statisticButton;
	JButton okButton;
	Color defaultColor;
	JFileChooser fc;
	HashMap<String, String> statisticmethod;
	public static Color buttonColor = new Color(255, 246, 143);

	JFrame frame;
	StatisticList StatisticList;
//	Vector<String> names = new Vector<String>();
	Color bColor;

	private JButton typeHButton = new JButton(Util.createImageIcon("Help16.gif"));

	/**
	 * Constructors the Repeat interface
	 */
	public ListDialog(JFrame frame, Vector<String> data,
			JButton mainButton, Color bColor, Color defaultColor,
			JFileChooser fc) {
		super(frame, "Repeat Data Files", false);

		this.data = data;
		this.statisticlist = new Vector<String>();
		this.statisticmethod = new HashMap<String, String>();
		this.frame = frame;
		this.mainButton = mainButton;
		this.bColor = bColor;
		this.defaultColor = defaultColor;
		this.fc = fc;
		okButton = new JButton("OK");
		okButton.setActionCommand("ok");
		okButton.addActionListener(this);

		okButton.setPreferredSize(new Dimension(145, 28));
		okButton.setMaximumSize(new Dimension(145, 28));
		removeButton = new JButton("Remove File");
		removeButton.setPreferredSize(new Dimension(145, 28));
		removeButton.setMaximumSize(new Dimension(145, 28));
		removeButton.setActionCommand("remove");
		viewButton = new JButton("View Selected File");
		viewButton.setPreferredSize(new Dimension(145, 28));
		viewButton.setMaximumSize(new Dimension(145, 28));
		viewButton.setActionCommand("view");
		statisticButton = new JButton("Statistic Method");
		statisticButton.setPreferredSize(new Dimension(145, 28));
		statisticButton.setMaximumSize(new Dimension(145, 28));
		statisticButton.setActionCommand("statistic");
		removeButton.addActionListener(this);
		removeButton.setEnabled((data.size() > 0));
		viewButton.addActionListener(this);
		viewButton.setEnabled((data.size() > 0));
		statisticButton.addActionListener(this);
		statisticButton.setEnabled((data.size() > 0));
		typeHButton.addActionListener(this);
		typeHButton.setActionCommand("help");

		final JButton setButton = new JButton("Add File",Util.createImageIcon("Open16.gif"));
		setButton.setActionCommand("add");
		setButton.addActionListener(this);
		setButton.setPreferredSize(new Dimension(145, 28));
		setButton.setMaximumSize(new Dimension(145, 28));
		getRootPane().setDefaultButton(setButton);

		// main part of the dialog
		list = new JList(data) {
			// Note from the sample
			// Subclass JList to workaround bug 4832765, which can cause the
			// scroll pane to not let the user easily scroll up to the beginning
			// of the list. An alternative would be to set the unitIncrement
			// of the JScrollBar to a fixed value. You wouldn't get the nice
			// aligned scrolling, but it should work.
			public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
				int row;
				if (orientation == SwingConstants.VERTICAL && direction < 0
						&& (row = getFirstVisibleIndex()) != -1) {
					Rectangle r = getCellBounds(row, row);
					if ((r.y == visibleRect.y) && (row != 0)) {
						Point loc = r.getLocation();
						loc.y--;
						int prevIndex = locationToIndex(loc);
						Rectangle prevR = getCellBounds(prevIndex, prevIndex);

						if (prevR == null || prevR.y >= r.y) {
							return 0;
						}
						return prevR.height;
					}
				}
				return super.getScrollableUnitIncrement(visibleRect, orientation, direction);
			}
		};
		list.setVisibleRowCount(-1);
		list.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					setButton.doClick(); // emulate button click
				}
			}
		});
		list.setSelectedIndex(0);
		JScrollPane listScroller = new JScrollPane(list);
		listScroller.setPreferredSize(new Dimension(300, 150));
		listScroller.setAlignmentX(LEFT_ALIGNMENT);

		
		list2 = new JList(statisticlist) {
			public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation, int direction) {
				int row;
				if (orientation == SwingConstants.VERTICAL && direction < 0
						&& (row = getFirstVisibleIndex()) != -1) {
					Rectangle r = getCellBounds(row, row);
					if ((r.y == visibleRect.y) && (row != 0)) {
						Point loc = r.getLocation();
						loc.y--;
						int prevIndex = locationToIndex(loc);
						Rectangle prevR = getCellBounds(prevIndex, prevIndex);
						if (prevR == null || prevR.y >= r.y) {
							return 0;
						}
						return prevR.height;
					}
				}
				return super.getScrollableUnitIncrement(visibleRect, orientation, direction);
			}
		};
		list2.setVisibleRowCount(-1);
		list2.addMouseListener(new MouseAdapter() {
			public void mouseClicked(MouseEvent e) {
				if (e.getClickCount() == 2) {
					setButton.doClick(); // emulate button click
				}
			}
		});
		list2.setSelectedIndex(0);
		JScrollPane listScroller2 = new JScrollPane(list2);
		listScroller2.setPreferredSize(new Dimension(300, 150));
		listScroller2.setAlignmentX(LEFT_ALIGNMENT);
		
		// Create a container so that we can add a title around
		// the scroll pane. Can't add a title directly to the
		// scroll pane because its background would be white.
		// Lay out the label and scroll pane from top to bottom.
		JPanel listPane = new JPanel();
		listPane.setLayout(new BoxLayout(listPane, BoxLayout.PAGE_AXIS));
		JLabel label = new JLabel("Repeat Data File(s):");
		label.setLabelFor(list);
		listPane.add(label);
		listPane.add(listScroller);
		listPane.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		listPane.setAlignmentY(TOP_ALIGNMENT);
		
		JPanel listPane2 = new JPanel();
		listPane2.setLayout(new BoxLayout(listPane2, BoxLayout.PAGE_AXIS));
		JLabel label2 = new JLabel("Statistic Method(s):");
		label2.setLabelFor(list2);
		listPane2.add(label2);
		listPane2.add(listScroller2);
		listPane2.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));
		listPane2.setAlignmentY(TOP_ALIGNMENT);
		
		JPanel listPane3 = new JPanel();
		listPane3.setLayout(new BoxLayout(listPane3, BoxLayout.X_AXIS));
		listPane3.add(listPane);
		listPane3.add(listPane2);
		listPane3.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));
		listPane3.setBackground(bColor);
		//listPane3.setAlignmentY(CENTER_ALIGNMENT);
		
		// Lay out the buttons from left to right.
		JPanel buttonPane = new JPanel(new SpringLayout());
		buttonPane.add(setButton);
		buttonPane.add(removeButton);
		buttonPane.add(viewButton);
		buttonPane.add(statisticButton);
		buttonPane.setBackground(bColor);
		SpringUtilities.makeCompactGrid(buttonPane, 1, 4, 6, 6, 6, 6);
		
		JPanel buttonPane2 = new JPanel(new SpringLayout());
		buttonPane2.add(okButton);
		buttonPane2.add(typeHButton);
		buttonPane2.setBackground(bColor);
		SpringUtilities.makeCompactGrid(buttonPane2, 1, 2, 6, 6, 6, 6);

		
		
		// Put everything together, using the content pane's BorderLayout.
		Container contentPane = getContentPane();
		BoxLayout layout = new BoxLayout(contentPane, BoxLayout.PAGE_AXIS);
		contentPane.setLayout(layout);
		contentPane.setBackground(bColor);
		// Initialize values.
		contentPane.add(listPane3);
		contentPane.add(buttonPane);
		contentPane.add(buttonPane2);

		pack();
	}

	/**
	 * Handle clicks on the Set and Cancel buttons.
	 */
	public void actionPerformed(ActionEvent e) {
		String szcommand = e.getActionCommand();
		if (szcommand.equals("ok")) {
			setVisible(false);
		} else if (szcommand.equals("view")) {
			Object[] nAselectindex = list.getSelectedValues();
			for (int aindex = 0; aindex < nAselectindex.length; aindex++) {
				final Object selectedval = nAselectindex[aindex];
				if (selectedval != null) {
					final JFrame fframe = frame;
					final String szfile = selectedval.toString();
					if ((new File(szfile)).exists()) {
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								DataTable newContentPane = new DataTable(fframe, szfile, false);
								JFrame dtframe = new JFrame(szfile);
								dtframe.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								dtframe.setLocation(20, 50);
								newContentPane.setOpaque(true); // content panes
																// must be
																// opaque
								dtframe.setContentPane(newContentPane);
								// Display the window.
								dtframe.pack();
								dtframe.setVisible(true);
							}
						});
					} else {
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, "File '"
										+ szfile + "' was not found.", "Error",
										JOptionPane.ERROR_MESSAGE);
							}
						});
					}
				}
			}
		} else if (szcommand.equals("statistic")) {
			statisticlist.clear();
			Object[] nAselectindex = list.getSelectedValues();
			for (int aindex = 0; aindex < nAselectindex.length; aindex++) {
				final Object selectedval = nAselectindex[aindex];
				if (selectedval != null) {
					final JFrame fframe = frame;
					final String szfile = selectedval.toString();
					if ((new File(szfile)).exists()) {

						
						statisticmethod.put(szfile, szfile+"11");
							
						StatisticList = new StatisticList(frame, this.bColor, this.defaultColor);
						StatisticList.setLocation(this.getX() + 75, this.getY() + 100);
						StatisticList.setVisible(true);
						
						
					} else {
						javax.swing.SwingUtilities.invokeLater(new Runnable() {
							public void run() {
								JOptionPane.showMessageDialog(fframe, "File '"
										+ szfile + "' was not found.", "Error",
										JOptionPane.ERROR_MESSAGE);
							}
						});
					}
				}
			}
			if(data.size()>0){
				for(int i=0;i<data.size();i++){
					if(statisticmethod.get(data.get(i)) != null){
						statisticlist.add(statisticmethod.get(data.get(i)));
					}else{
						statisticlist.add(" ");
					}
				}
			}
			
		} else if (szcommand.equals("remove")) {
			Object[] nAselectindex = list.getSelectedValues();
			for (int aindex = 0; aindex < nAselectindex.length; aindex++) {
				Object selectedval = nAselectindex[aindex];
				if (selectedval != null) {
					data.remove(selectedval);
					statisticlist.remove(statisticmethod.get(selectedval));
				}
			}

			if (data.size() < 1) {
				removeButton.setEnabled(false);
				viewButton.setEnabled(false);
				statisticButton.setEnabled(false);
				mainButton.setBackground(defaultColor);
			} else if (list.getSelectedIndex() >= data.size()) {
				list.setSelectedIndex(data.size() - 1);
			}
		} else if (szcommand.equals("add")) {
			int returnVal = fc.showOpenDialog(this);

			if (returnVal == JFileChooser.APPROVE_OPTION) {
				File file = fc.getSelectedFile();
				data.add(file.getAbsolutePath());
				removeButton.setEnabled(true);
				viewButton.setEnabled(true);
				statisticButton.setEnabled(true);
				mainButton.setBackground(ListDialog.buttonColor);
				list.setSelectedIndex(data.size() - 1);
			}
		} else if (szcommand.equals("help")) {
			String szMessage = " Click 'Add File' to select a repeat file to add "
					+ "to the list of repeat files.  To remove a file click on a file from the list "
					+ "and then click 'Remove File'.  To view a file click on the file and then "
					+ "click 'View Selected File'.  Click 'OK' to close "
					+ "the dialog window.\n\n"
					+ "Repeat files "
					+ "need to be in the same format as the original data file with "
					+ "the same genes in the same order.  "
					+ "All gene expression values will be normalized against their expression at first time "
					+ "point will be computed for each data set (repeat(s) and original).\n\n"
					+ "Repeat data can either represent repeat measurements taken on the same samples as the "
					+ "original experiment ('The Same Time Period'), or distinct full repeat experiments from "
					+ "different time periods ('Different Time Periods').  "
					+ "If the repeat data is from different time periods genes will be filtered that do not display "
					+ "a consistent expression profile between repeats "
					+ "(based on a parameter under the advanced options).  "
					+ "Also if the repeat data is from different time periods, then the clustering of model profiles "
					+ "can be adjusted using noise estimates from the repeat experiments (also based on a parameter "
					+ "under the advanced options).\n\n"
					+ "If the values are from 'The Same Time Period' then each time point value is averaged by the "
					+ "median (in the case of two time points a simple mean is used).  If the values are from "
					+ "'Different Time Periods' the median of values after normalization is used.\n\n";

			Util.renderDialog(frame, szMessage);
		}
		list.updateUI();
		list2.updateUI();
	}

	/**
	 * Updates the repeat list interface based on the contents of updateData and
	 * bupdatealltime
	 */
	public void updateSettings(Vector<String> updateData, boolean bupdatealltime) {

		list.setSelectedIndex(0);
		if (updateData.size() > 0) {
			removeButton.setEnabled(true);
			viewButton.setEnabled(true);
			mainButton.setBackground(ListDialog.buttonColor);
		} else {
			removeButton.setEnabled(false);
			viewButton.setEnabled(false);
			mainButton.setBackground(defaultColor);
		}

		data.removeAllElements();
		data.addAll(updateData);
		list.updateUI();
	}
}
