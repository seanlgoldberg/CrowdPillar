package iitb.CRF;

import java.util.ArrayList;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;

public class Compute {
	
	ArrayList<double[]> Alpha_i;
	double[][] Mi;
	int seqlen;
	DoubleMatrix1D beta[];
	boolean initMDone;
	
	public Compute(int seqlength) {
		
		Alpha_i=new ArrayList();
		double alpha[]=new double[7];
		alpha[0]=1;
		for (int i = 1; i <7; i++)
			alpha[i] = 0;
		Alpha_i.add(alpha);
		seqlen=seqlength;
	}
	
	public void initializeMi(DoubleMatrix2D Mi) {
		this.Mi=new double[7][7]; 
		for (int i = 0; i < Mi.rows(); i++)
			for (int j = 0; j < Mi.columns(); j++)
				this.Mi[i][j] = Mi.get(i, j);
	}
	
	public void updateAlphaOne(DoubleMatrix1D Ri) {
		double alpha[]=new double[7];
		for (int i = 0; i < Ri.size(); i++)
			alpha[i]=Ri.get(i);
		Alpha_i.add(alpha);
	}
	
	public void generateAlphaAtNode(int nodeNo) {
		double newAlpha[] = new double[7];
		double oldAlpha[] = (double[]) Alpha_i.get(nodeNo-1);
		for (int i = 0; i < 7; i++) {
			for (int j = 0; j < 7; j++) {
				newAlpha[i] += oldAlpha[j] * Mi[j][i];
			}
			
		}
		Alpha_i.add(newAlpha);
	}
	
	public void printAlpha(int i) {
		double[] d= (double[]) Alpha_i.get(i);
		System.out.print("Alpha "+i+": ");
		for(int i1=0;i1<7;i1++)
			System.out.print(d[i1]+" ");
		System.out.println("\n");
	}
	
	
	

}
