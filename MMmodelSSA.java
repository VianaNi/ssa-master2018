
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
 *This model is specifically created for the Michalies-menten kinetics. 
 *Users can change fundamental element values in the code.
 *While only maximum simulation time initial states can be defined through the console.
 */
public class MMmodelSSA extends JFrame {

	private static final long serialVersionUID = 1L;
    static double matrixPre[][];
    static double matrixPost[][];
    int m;// the number of reactions
    int n;// the number of species
	static double[] h;// the hazard for each reaction
	static double[] rate;
	static int[] state;
	static double Tmax;// the overall reaction time
    static String[] strs;
	XYSeriesCollection dataset;


	public MMmodelSSA(String s) {
		super(s);
		this.m =3;
		this.n =4;
	}
	

	static Scanner sc = new Scanner(System.in);

	public static void main(String args[]) throws Exception {


        matrixPre  = new double [3][4];
		matrixPre [0][0]= 1;
		matrixPre [1][0]= 0;
		matrixPre [2][0]= 0;
		matrixPre [0][1]= 1;
		matrixPre [1][1]= 0;
		matrixPre [2][1]= 0;
		matrixPre [0][2]= 0;
		matrixPre [1][2]= 1;
		matrixPre [2][2]= 1;
		matrixPre [0][3]= 0;
		matrixPre [1][3]= 0;
		matrixPre [2][3]= 0;
		
	    matrixPost  = new double [3][4];
	    matrixPost [0][0]= 0;
		matrixPost [1][0]= 1;
		matrixPost [2][0]= 0;
		matrixPost [0][1]= 0;
		matrixPost [1][1]= 1;
		matrixPost [2][1]= 1;
		matrixPost [0][2]= 1;
		matrixPost [1][2]= 0;
		matrixPost [2][2]= 0;
		matrixPost [0][3]= 0;
		matrixPost [1][3]= 0;
		matrixPost [2][3]= 1;
		
		strs = new String [4];
		strs[0]="S";
		strs[1]="E";
		strs[2]="SE";
		strs[3]="P";
		
		rate = new double [3];
		rate[0]=0.00166;
		rate[1]=0.001;
		rate[2]=0.1;

		RealMatrix PreMatrix = new Array2DRowRealMatrix(matrixPre);

		RealMatrix PostMatrix = new Array2DRowRealMatrix(matrixPost);

		//output the stochiometry matrix
		RealMatrix matrixSubtract = PostMatrix.subtract(PreMatrix);

		// input the default maxium time of the whole reaction
		System.out.println("Enter the maxium time of the whole reaction:");
		double Tmax = sc.nextDouble();
		
		//
		double newTime =0;
		Vector<Double> timeList = new Vector<Double>(0);
		Vector<int[]> allStateList = new Vector<int[]>(0);
		timeList.add(newTime);

		// input the initial states of numbers of molecules
		System.out.println("Enter the initial molecule numbers of all species:");
		int[] state = new int[4];
		for (int i = 0; i < 4; i++) {
			state[i] = sc.nextInt();
		}
		allStateList.add(state);
		
		long startTime = System.currentTimeMillis(); 
	
		MARKONE:
		while (newTime < Tmax) {
			for (int loopIndex = 0; loopIndex < allStateList.size(); loopIndex++) {
				// calculate the hazard for each reaction and the general hazard
	                double [] bc =new double[3];
	                Arrays.fill(bc,1);
	                double [] h = new double [3];
	                 for (int i = 0; i < 3; i++) {
						for (int j = 0; j < 4; j++) {
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
				for (int i =0; i<3; i++) {
						H = DoubleStream.of(h).sum();
				}
				
				if(H==0) {
					break;
					}

				// select a random reaction time
				double tau = -(Math.log(Math.random()))/H;


				// select a random reaction j with discrete random probability;
				int[] indexList = new int[3];
				double [] probabilities = new double[3];
		
               for (int i =0; i<3; i++) {
					indexList[i]=i;
					probabilities[i]=h[i]/H;
				}   
          EnumeratedIntegerDistribution EID = new EnumeratedIntegerDistribution(indexList,probabilities);
				int index = EID.sample();

				// Update the state
				double[] vectorDoubleArray = matrixSubtract.transpose().getColumnVector(index).toArray();
				int[] intArray = new int[4];
				for (int i = 0; i < 4; i++) {
					intArray[i] = (int) vectorDoubleArray[i];
				}
				int[] newState = new int[4];
				for (int p = 0; p < 4; p++) {
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
		long endTime   = System.currentTimeMillis();
		long TotalTime = endTime - startTime; 
		System.out.println("So the computer used "+ TotalTime +" milliseconds to compute the SSA model");
		
	
		// close the scanner
		sc.close();

       // form the dataset
		XYSeriesCollection dataset = new XYSeriesCollection();
		XYSeries[] species = new XYSeries[4];
			species[0]=new XYSeries("S");
			species[1]=new XYSeries("E");
			species[2]=new XYSeries("SE");
			species[3]=new XYSeries("P");
	
		for (int j =0; j<4; j++) {
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
                "SSA Modeling of M-M Kinetics", 
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
	        for (int j =0; j<4;j++){
	        
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

}

