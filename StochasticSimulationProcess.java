
import java.util.*;

import java.util.stream.DoubleStream;
import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Font;
import java.lang.Math;


import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.block.BlockBorder;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.TextTitle;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


/**
 * @author zhang nijia
 *This prototype is created especially for the stochastic simulation algorithm.
 */
public class StochasticSimulationProcess extends JFrame {
	
	private static final long serialVersionUID = 1L;
	static int matrixPre[][];
	static int matrixPost[][];
	static int m;// the number of reactions
	static int n;// the number of species
	static double[] rate;
	static int[] state;
	static double Tmax;// the overall reaction time
	static String[] strs;
	XYSeriesCollection dataset;


	public StochasticSimulationProcess(String s) {
		super(s);
	}

	static Scanner sc = new Scanner(System.in);

	public static void main(String args[]) throws Exception {

		// input the number of species
		System.out.println("Enter the number of species:");
		int n = sc.nextInt();

		// input the number of reactions
		System.out.println("Enter the number of reactions:");
		int m = sc.nextInt();

		//enter the pre and post part of coefficients
		int[][] matrixPre = new int[m][n];
		enterMatrixDataPre(sc, matrixPre, m, n);
		printMatrixPre(matrixPre, m, n);

		// convert the 2d int to 2d double
		double[][] matrixPre2 = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				matrixPre2[i][j] = (double) matrixPre[i][j];
		}

		RealMatrix PreMatrix = new Array2DRowRealMatrix(matrixPre2);

		// remember to add space key when doing the typing
		int[][] matrixPost = new int[m][n];
		enterMatrixDataPost(sc, matrixPost, m, n);
		printMatrixPost(matrixPost, m, n);

		// convert the 2d int to 2d double
		double[][] matrixPost2 = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++)
				matrixPost2[i][j] = (double) matrixPost[i][j];
		}

		RealMatrix PostMatrix = new Array2DRowRealMatrix(matrixPost2);

		//output the stochiometry matrix
		RealMatrix matrixSubtract = PostMatrix.subtract(PreMatrix);

		// input the default maxium time of the whole reaction
		System.out.println("Enter the maxium simulation time:");
		double Tmax = sc.nextDouble();

		// input the name of all the species
		System.out.println("Enter the name of all species");
		String[] strs = new String[n];
		for (int i = 0; i < n; i++) {
			strs[i] = sc.nextLine();
		}

		//input the rate constant
		System.out.println("Enter the rate constant of each reaction:");
		Double[] rate = new Double[m];
		for (int r = 0; r < m; r++) {
			rate[r] = sc.nextDouble();
		}

		//
		double newTime =0;
		Vector<Double> timeList = new Vector<Double>(0);
		Vector<int[]> allStateList = new Vector<int[]>(0);
		timeList.add(newTime);

		// input the initial states of numbers of molecules
		System.out.println("Enter the initial molecule numbers of all species:");
		int[] state = new int[n];
		for (int i = 0; i < n; i++) {
			state[i] = sc.nextInt();
		}
		allStateList.add(state);
		
		long startTime = System.currentTimeMillis(); 
		MARKONE:
		while (newTime < Tmax) {
			for (int loopIndex = 0; loopIndex < allStateList.size(); loopIndex++) {
				// calculate the propensity for each reaction and the general propensity
	                double [] bc =new double[m];
	                Arrays.fill(bc,1);
	                double [] h = new double [m];
	                 for (int i = 0; i < m; i++) {
						for (int j = 0; j < n; j++) {
							if (allStateList.get(loopIndex)[j]< (int) (PreMatrix.getRowVector(i).toArray()[j])) {
								bc[i]=0;	
							}	
							else {
							bc[i]*=CombinatoricsUtils.binomialCoefficient(
								allStateList.get(loopIndex)[j], (int) (PreMatrix.getRowVector(i).toArray()[j]));
							}	
						}
						 h[i] = rate[i]*bc[i];
						
	                 }
								
				double H =0;
				for (int i =0; i<m; i++) {
						H = DoubleStream.of(h).sum();
				}
				
				if(H==0) {
					break;
					}

				// select a random reaction time
				double tau = -(Math.log(Math.random()))/H;

				// select a random reaction j with discrete random probability;
				int[] indexList = new int[m];
				double [] probabilities = new double[m];
		
               for (int i =0; i<m; i++) {
					indexList[i]=i;
					probabilities[i]=h[i]/H;
				}   
          EnumeratedIntegerDistribution EID = new EnumeratedIntegerDistribution(indexList,probabilities);
				int index = EID.sample();
				System.out.println("So the next simulated event is:" + index);

				// Update the state
				double[] vectorDoubleArray = matrixSubtract.transpose().getColumnVector(index).toArray();
				int[] intArray = new int[n];
				for (int i = 0; i < n; i++) {
					intArray[i] = (int) vectorDoubleArray[i];
				}
				int[] newState = new int[n];
				for (int p = 0; p < n; p++) {
					newState[p] = intArray[p] + allStateList.get(loopIndex)[p];
				}
				allStateList.add(newState);
				
				// put the newTime
				newTime += tau;
				timeList.add(newTime);
				
				System.out.printf("%-3s %s\n", loopIndex, newTime+ "   " + Arrays.toString(newState));
				
				if(newTime>=Tmax) {
					System.out.println("So the simulation time is reached") ;
					break MARKONE;
				}
			}
			
		}
		long endTime = System.currentTimeMillis();
		long TotalTime = endTime - startTime; 
		System.out.println("So the computer used "+ TotalTime +" milliseconds to compute the SSA model");
		
	
		// close the scanner
		sc.close();

       // form the dataset
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries[] species = new XYSeries[n];
		for (int j =0; j<n; j++) {
			species[j]=new XYSeries(""+Integer.toString(j));
		}
		for (int j =0; j<n; j++) {
			for (int k=0; k< allStateList.size();k++) {
			species[j].add(timeList.get(k).doubleValue(), allStateList.get(k)[j]);
		}
			dataset.addSeries(species[j]);
	}
           
		// plot out the graph

		boolean showLegend = true;
		boolean createURL = true;
		boolean createTooltip = true;
		
		JFreeChart chart = ChartFactory.createXYLineChart(
                "SSA Modeling", 
                "time", 
                "Number of Molecules", 
                dataset, 
                PlotOrientation.VERTICAL,
                showLegend, 
                createURL, 
                createTooltip
        );

		 ChartPanel chartPanel = new ChartPanel(chart);
	        chartPanel.setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));
	        chartPanel.setBackground(Color.white);
	       XYPlot plot = chart.getXYPlot();
	     
	       
	        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
	        for (int j =0; j<n;j++){
	        
	       // get random colors for each XYSeries
	     renderer.setSeriesPaint(j, Color.getHSBColor(0 + (float)(Math.random() * ((255 - 0) + 1)), 0 + (float)(Math.random() * ((255 - 0) + 1)), 0 + (float)(Math.random() * ((255 - 0) + 1))));
	         renderer.setSeriesStroke(j, new BasicStroke(2.0f));
	        }
	      
	     
            renderer.setBaseItemLabelsVisible(true);

	        plot.setRenderer(renderer);
	        plot.setBackgroundPaint(Color.white);

	        plot.setRangeGridlinesVisible(false);
	        plot.setDomainGridlinesVisible(false);
             
	       
	       chart.getLegend().setFrame(BlockBorder.NONE);

	        chart.setTitle(new TextTitle("Stochastic Simulation Algorithm",
	                        new Font("Serif", Font.BOLD, 18)
	                )
	        );

		SwingUtilities.invokeLater(() -> {
			StochasticSimulationProcess ex = new StochasticSimulationProcess("SSA");
	        ex.setContentPane(chartPanel);
			ex.pack();
	        ex.setLocationRelativeTo(null);
	        ex.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);  
			ex.setVisible(true);
	
		});
		
	}

	public static void enterMatrixDataPre(Scanner sc, int[][] matrixPre, int m, int n) {
		System.out.println("Enter Matrix pre data");
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				matrixPre[i][j] = sc.nextInt();
			}
		}
	}

	public static void printMatrixPre(int[][] matrixPre, int m, int n) {
		System.out.println("Matrix pre is:");
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				System.out.print(matrixPre[i][j] + "\t");
			}
			System.out.println();
		}
	}

	public static void enterMatrixDataPost(Scanner sc, int[][] matrixPost, int m, int n) {
		System.out.println("Enter Matrix post data");
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				matrixPost[i][j] = sc.nextInt();
			}
		}
	}

	public static void printMatrixPost(int[][] matrixPost, int m, int n) {
		System.out.println("Matrix post is:");
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				System.out.print(matrixPost[i][j] + "\t");
			}
			System.out.println();
		}
	}
}
