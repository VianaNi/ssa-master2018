
import java.awt.BasicStroke;

import java.awt.Color;
import java.awt.Font;
import java.util.Arrays;
import java.util.Scanner;
import java.util.Vector;
import java.util.stream.DoubleStream;

import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.SwingUtilities;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.distribution.EnumeratedIntegerDistribution;
import org.apache.commons.math3.distribution.PoissonDistribution;
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
 *The model is created especially for the poisson tau-leap algroithm
 */
public class PoissonTauLeap extends JFrame{
	private static final long serialVersionUID = 1L;
	static int matrixPre[][];
	static int matrixPost[][];
	static int m;// the number of reactions
	static int n;// the number of species
	static double[] h;// the hazard for each reaction
	static double[] rate;
	static int max;
	static double Tmax;// the overall reaction time
	Vector<Double> timeList;
	Vector<int[]> allStateList;
	static String[] strs;
	XYSeriesCollection dataset;
	
	public PoissonTauLeap(String s) {
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

		//
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

		//
		RealMatrix matrixSubtract = PostMatrix.subtract(PreMatrix);
		System.out.println("So the transpose matrix after subtraction is:\t" + matrixSubtract.transpose());

		// input the default maxium time of the whole reaction
		System.out.println("Enter the maxium time of the whole reaction:");
		double Tmax = sc.nextDouble();

		// input the name of all the species
		System.out.println("Enter the name of all species");
		String[] strs = new String[n];
		for (int i = 0; i < n; i++) {
			strs[i] = sc.nextLine();
		}

		// input the rate constant of all the reactions(the coefficient), must be a
		// double
		System.out.println("Enter the rate constant of each reaction:");
		Double[] rate = new Double[m];
		for (int r = 0; r < m; r++) {
			rate[r] = sc.nextDouble();
		}

		//
		Vector<Double> timeList = new Vector<Double>(0);
		Vector<int[]> allStateList = new Vector<int[]>(0);
		double newTime =0;
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
  while(newTime<Tmax) {
	for (int loopIndex = 0; loopIndex < allStateList.size(); loopIndex++) {
		// calculate the hazard for each reaction and the general hazard
        double [] bc =new double[m];
        Arrays.fill(bc,1);
        double [] h = new double [m];
         for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (allStateList.get(loopIndex)[j]<(int) (PreMatrix.getRowVector(i).toArray()[j])) {
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
		break MARKONE;
		}

	
		
	// verify the -th order reaction and define gi
	int [] g = new int[n];
	for (int i =0; i<n; i++) {
		for (int j =0; j<m; j++) {
			 double max = PreMatrix.getColumnVector(i).toArray()[0];
			if(PreMatrix.getColumnVector(i).toArray()[j]>max) {
				max = PreMatrix.getColumnVector(i).toArray()[j];
			}
			if (max ==1) {
				g[i] = 1;
			}
			else if (max ==2){
				if(PreMatrix.getColumnVector(i).toArray()[j]==2) {
				try{
					g[i] = 2+1/(allStateList.get(loopIndex)[i]-1);
					}
				catch (java.lang.ArithmeticException exception){
					g[i]=1;
				}
				}
				else if (PreMatrix.getColumnVector(i).toArray()[j]!=2) {
				 g[i] = 2;
				}
			}
			else if (max == 3) {
				if(PreMatrix.getColumnVector(i).toArray()[j]==2) {
					g[i] = (2+1/(allStateList.get(loopIndex)[i]-1))*(3/2);
					}
				else if (PreMatrix.getColumnVector(i).toArray()[j]==3) {
					 g[i] = 3+1/(allStateList.get(loopIndex)[i]-1)+2/(allStateList.get(loopIndex)[i]-2);
			}
				else if(PreMatrix.getColumnVector(i).toArray()[j]!=3 && PreMatrix.getColumnVector(i).toArray()[j]!=2) {
					g[i] = 3;
				}
		   }
			else if (max == 0) {
				g[i] = 0;
			}
		}
	}
	
	// simulate the time step
	// Cao defined that sigma equals to 0.03.
	
	double expectedValue[]= new double[n];
	double variance[] = new double[n];
	double tauList[] = new double[n];
	double tau = tauList[0];
	
	    for(int i =0; i<n; i++) {
	    	for (int j =0; j<m; j++) {
	expectedValue[i]+= h[j]*(matrixSubtract.getRowVector(j).toArray()[i]);
	variance[i] += h[j]*FastMath.pow((matrixSubtract.getRowVector(j).toArray()[i]),2);
	    }
	}
	    for (int i =0; i<n; i++) {
	    	if(g[i]==0) {
	        tauList[i]=0.01;
	        }
	        else {
	        	tauList[i]= Math.min((Math.max((0.03/g[i])*allStateList.get(loopIndex)[i], 1)/FastMath.abs(expectedValue[i])),
	        			(FastMath.pow(Math.max((0.03/g[i])*allStateList.get(loopIndex)[i],1), 2)/variance[i]));
	        }
	    	if(tau == 0 || tauList[i]<tau)
	    		tau = tauList[i];	
	    }
	 newTime +=tau; 

	
	//affirm the firing times
	int p = 1;
	int [] kList = new int [m];
    PoissonDistribution poisson = new PoissonDistribution(p);
    int [] change = new int [n];
    int [] newState = new int[n];
    
    for (int j =0;j<m; j++) { 
   	 p =(int)(h[j]*tau);
	    int k = (int)poisson.sample();	
	    kList[j] =k;
    }
    System.out.println("So the firing times are:"+Arrays.toString(kList));
    
    
 
    for (int i =0; i<n; i++) {
	    for (int j =0;j<m; j++) { 
    		double[] vectorDoubleArray = matrixSubtract.getColumnVector(i).toArray();
    	    change [i]+= kList[j]*(int) vectorDoubleArray[j];
    	    newState[i]= allStateList.get(loopIndex)[i]+change[i];
	    }
    }
 
    for (int i =0; i<n; i++) {
    		while(newState[i]<0) {
    		newTime = newTime-tau;
	    	  tau = -(Math.log(Math.random()))/H;	  
	          int[] indexList = new int[m];
				double [] probabilities = new double[m];
             for (int b =0; b<m; b++) {
					indexList[b]=b;
					probabilities[b]=h[b]/H;
				}     
        EnumeratedIntegerDistribution EID = new EnumeratedIntegerDistribution(indexList,probabilities);
				int index = EID.sample(); 
				 double[] vectorDoubleArray = matrixSubtract.transpose().getColumnVector(index).toArray();
				 for (int b=0; b<n; b++) {
					 if(loopIndex>0) {
				newState[b] =(int) vectorDoubleArray[b] + allStateList.get(loopIndex-1)[b];
					 }   
					 else if(loopIndex==0) {
				 newState[b] =(int) vectorDoubleArray[b] +state[b]; 
				 }
	       } 
				 newTime+=tau;
       }
    
  
	if(newTime>=Tmax) {
		System.out.println("So the simulation time is reached") ;
		break MARKONE;
	   }
      } 
    allStateList.add(newState);
    timeList.add(newTime);
    System.out.printf("%-3s %s\n", loopIndex, newTime+ "   " + Arrays.toString(newState));
	}
  }
		long endTime   = System.currentTimeMillis();
		long TotalTime = endTime - startTime; 
		System.out.println("So the computer used "+ TotalTime + " milliseconds to compute the tau-leap model");
	// close the scanner
			sc.close();

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
		
		JFreeChart chart = ChartFactory.createXYLineChart(
                "Poisson Tau Leap Algorithm", 
                "time", 
                "Number of Molecules", 
                dataset, 
                PlotOrientation.VERTICAL,
                true, 
                true, 
                false
        );

		 ChartPanel chartPanel = new ChartPanel(chart);
	        chartPanel.setBorder(BorderFactory.createEmptyBorder(15, 15, 15, 15));
	        chartPanel.setBackground(Color.white);
	       XYPlot plot = chart.getXYPlot();
	     
	       

	        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer();
	        for (int j =0; j<n;j++){
	        
	       // get random colors for each XYSeries
	         renderer.setSeriesPaint(j, Color.getHSBColor(0 + (float)(Math.random() * ((255 - 0) + 1)), 
	        		 0 + (float)(Math.random() * ((255 - 0) + 1)), 0 + (float)(Math.random() * ((255 - 0) + 1))));
	         renderer.setSeriesStroke(j, new BasicStroke(2.0f));
	        }
	               
	        plot.setRenderer(renderer);
	        plot.setBackgroundPaint(Color.white);

	        plot.setRangeGridlinesVisible(false);
	        plot.setDomainGridlinesVisible(false);

	        chart.getLegend().setFrame(BlockBorder.NONE);

	        chart.setTitle(new TextTitle("Poisson Tau Leap Algorithm",
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
}}

