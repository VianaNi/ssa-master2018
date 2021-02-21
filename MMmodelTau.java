
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
 *This model is specifically created for the Michalies-menten kinetics. 
 *Users can change fundamental element values in the code.
 *While only maximum simulation time initial states can be defined through the console.
 */
public class MMmodelTau extends JFrame {
	private static final long serialVersionUID = 1L;
	static double matrixPre[][];
	static double matrixPost[][];
    int m;// the number of reactions
    int n;// the number of species
	static double[] h;// the hazard for each reaction
	static double[] rate;
	static double Tmax;// the overall reaction time
	Vector<Double> timeList;
	Vector<int[]> allStateList;
	static String[] strs;
	XYSeriesCollection dataset;
	
	public MMmodelTau(String s) {
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
			rate[1]=0.00001;
			rate[2]=0.000005;

		RealMatrix PreMatrix = new Array2DRowRealMatrix(matrixPre);

		RealMatrix PostMatrix = new Array2DRowRealMatrix(matrixPost);

		//
		RealMatrix matrixSubtract = PostMatrix.subtract(PreMatrix);

		// input the default maxium time of the whole reaction
		System.out.println("Enter the maxium time of the whole reaction:");
		double Tmax = sc.nextDouble();


		//
		Vector<Double> timeList = new Vector<Double>(0);
		Vector<int[]> allStateList = new Vector<int[]>(0);
		double newTime =0;
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
  while(newTime<Tmax) {
	for (int loopIndex = 0; loopIndex < allStateList.size(); loopIndex++) {
		// calculate the hazard for each reaction and the general hazard
        double [] bc =new double[3];
        Arrays.fill(bc,1);
        double [] h = new double [3];
         for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 4; j++) {
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
	for (int i =0; i<3; i++) {
			H = DoubleStream.of(h).sum();
	}
	if(H==0) {
		break MARKONE;
		}

	double g[] = {1,1,1,0};

	
	// simulate the time step
	// Cao defined that sigma equals to 0.03.
	double expectedValue[]= new double[4];
	double variance[] = new double[4];
	 double tauList[] = new double[3];
	 double tau = tauList[0];
	 
	    for(int i =0; i<4; i++) {
	    	for (int j =0; j<3; j++) {
	expectedValue[i]+= h[j]*(matrixSubtract.getRowVector(j).toArray()[i]);
	variance[i] += h[j]*FastMath.pow((matrixSubtract.getRowVector(j).toArray()[i]),2);
	    }
	}
	    for (int i =0; i<3; i++) {
	        	tauList[i]= Math.min((Math.max((0.03/g[i])*allStateList.get(loopIndex)[i], 1)/FastMath.abs(expectedValue[i])),(FastMath.pow(Math.max((0.03/g[i])*allStateList.get(loopIndex)[i],1), 2)/variance[i]));
	        	if(tau == 0 || tauList[i]<tau)
		    		tau = tauList[i];	
	    }
	   
	 newTime +=tau; 

	
	//affirm the firing times
	int p = 1;
	int [] kList = new int [3];
    PoissonDistribution poisson = new PoissonDistribution(p);
    int [] change = new int [4];
    int [] newState = new int[4];
    
    for (int j =0;j<3; j++) { 
   	 p =(int)(h[j]*tau);
	    int k = (int)poisson.sample();	 
	    kList[j] =k;
    }
    
    for (int i =0; i<4; i++) {
	    for (int j =0;j<3; j++) { 
    		double[] vectorDoubleArray = matrixSubtract.getColumnVector(i).toArray();
    	    change [i]+= kList[j]*(int) vectorDoubleArray[j];
    	    newState[i]= allStateList.get(loopIndex)[i]+change[i];
	    }
    }
 
    for (int i =0; i<4; i++){
    		while(newState[i]<0) {
    		newTime = newTime-tau;
	    	  tau = -(Math.log(Math.random()))/H;	  
	          int[] indexList = new int[3];
				double [] probabilities = new double[3];
             for (int b =0; b<3; b++) {
					indexList[b]=b;
					probabilities[b]=h[b]/H;
				}     
        EnumeratedIntegerDistribution EID = new EnumeratedIntegerDistribution(indexList,probabilities);
				int index = EID.sample();
				 double[] vectorDoubleArray = matrixSubtract.transpose().getColumnVector(index).toArray();
				 for (int b=0; b<4; b++) {
					 if(loopIndex>0) {
				newState[b] =(int) vectorDoubleArray[b] + allStateList.get(loopIndex-1)[b];
					 }   
					 else if(loopIndex==0) {
				 newState[b] =(int) vectorDoubleArray[b] +state[b];
				 }
	       } 
       }
    		newTime+=tau;
    }
   
    allStateList.add(newState);
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
		System.out.println("So the computer used "+ TotalTime + " milliseconds to compute the tau-leap model");
	// close the scanner
			sc.close();

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
		
		JFreeChart chart = ChartFactory.createXYLineChart(
                "m-m kinetics Poisson Tau Leap Algorithm", 
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
	        for (int j =0; j<4;j++){
	        
	       // get random colors for each XYSeries
	         renderer.setSeriesPaint(j, Color.getHSBColor(0 + (float)(Math.random() * ((255 - 0) + 1)), 0 + (float)(Math.random() * ((255 - 0) + 1)), 0 + (float)(Math.random() * ((255 - 0) + 1))));
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

}


