
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import java.util.Scanner;
import java.util.Vector;
import java.util.List;


import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.util.CombinatoricsUtils;
import org.apache.commons.math3.util.FastMath;



/**
 * @author zhang nijia
 *This part of the code only attempts to offer a modification method to classify the critical
 *and non-critical reaction metnioned in the section7.4 in the paper. 
 *This part of the code doesn't represent the whole algorithm!!!!
 */
public class ModificationCode {
	static int matrixPre[][];
	static int matrixPost[][];
	static int m;// the number of reactions
	static int n;// the number of species
	static double[] h;// the hazard for each reaction
	static double[] rate;
	static int[] state;
	static double newTime;
	static int max;
	static double Tmax;// the overall reaction time
	Vector<Double> timeList;
	static double H;
	Vector<int[]> allStateList;
	static String[] strs;
	static double tau;
	
	
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
		RealMatrix matrixSubtract = PreMatrix.subtract(PostMatrix);
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
		timeList.add(newTime);

		// input the initial states of numbers of molecules
		System.out.println("Enter the initial molecule numbers of all species:");
		int[] state = new int[n];
		for (int i = 0; i < n; i++) {
			state[i] = sc.nextInt();
		}

		allStateList.add(state);

  while(newTime<Tmax) {
	for (int loopIndex = 0; loopIndex < allStateList.size(); loopIndex++) {
		// calculate the hazard
		double[] h = new double[m];
		H = 0;
		try {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					h[i] = rate[i] * CombinatoricsUtils.binomialCoefficientDouble(
							allStateList.get(loopIndex)[j], (int) (PreMatrix.getRowVector(i).toArray()[j]));
					H += h[i];
				}
			}
		}

		catch (NumberIsTooLargeException exceptionn) {
			System.out.println("One of the species has been exhausted and there is no next firing");
		}
		System.out.println("So the general hazard is:" + H);
		
	//
		// classify the critical and non-critical reaction
		List<Integer> criticalReactionList = new ArrayList<>();
		List<Integer> noncriReactionList = new ArrayList<>();
		double [] allSmallest = new double[m];	
		double [] result = new double[n];
	for (int j=0; j<m; j++) {
			for (int i =0; i<n; i++) {
				if(matrixSubtract.transpose().getColumnVector(j).toArray()[i]!=0) {
				result[i] = Math.floor(allStateList.get(loopIndex)[i]/FastMath.abs((int)matrixSubtract.transpose().getColumnVector(j).toArray()[i]));
				}
			// the denomiator cannot be zero			
				double smallestCR = result[0]; 
				if(result[i]<=smallestCR)
					allSmallest[j]=result[i];
			}
	
			if(allSmallest[j]>=10) {
			   criticalReactionList.add(j);
		  }
			else
				 noncriReactionList.add(j);
	}
	 
	Collections.sort(criticalReactionList);
	Collections.sort(noncriReactionList);
	
	//convert list to array
	int[] criticalReaction = new int[criticalReactionList.size()];
	for(int i = 0; i < criticalReactionList.size(); i++) criticalReaction[i] = criticalReactionList.get(i);
	
	int[] noncriReaction = new int[noncriReactionList.size()];
	for(int i = 0; i < noncriReactionList.size(); i++) noncriReaction[i] = noncriReactionList.get(i);
	
	
	System.out.println("So the critical reaction index are :" + Arrays.toString(criticalReaction));
	System.out.println("So the non-critical reaction index are :" + Arrays.toString(noncriReaction));
	
	//
	     }
     }
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

	