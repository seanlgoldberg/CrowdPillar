package iitb.CRF;

import iitb.CRF.Trainer.MultFunc;
import iitb.CRF.Trainer.MultSingle;

import java.io.Serializable;
import java.util.Vector;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.*;

/**
 * 
 * Viterbi search
 * 
 * @author Sunita Sarawagi
 * 
 */

public class Viterbi implements Serializable {
	private static final long serialVersionUID = 8122L;
	protected CRF model;
	protected int beamsize;
	protected boolean reuseM;
	EdgeGenerator edgeGen;

	Viterbi(CRF model, int bs) {
		this.model = model;
		beamsize = bs;
		if (model.params.miscOptions.getProperty("beamSize") != null)
			beamsize = Integer.parseInt(model.params.miscOptions
					.getProperty("beamSize"));

	}

	protected class Entry {
		public Soln solns[]; // TODO.
		boolean valid = true;

		protected Entry() {
		}

		protected Entry(int beamsize, int id, int pos) {
			solns = new Soln[beamsize];
			for (int i = 0; i < solns.length; i++)
				solns[i] = newSoln(id, pos);
			System.out.println(beamsize + " " + solns.length);
		}

		protected Soln newSoln(int label, int pos) {
			return new Soln(label, pos);
		}

		protected void clear() {
			valid = false;
			for (int i = 0; i < solns.length; i++)
				solns[i].clear();
		}

		protected int size() {
			return solns.length;
		}

		protected Soln get(int i) {
			return solns[i];
		}

		protected void insert(int i, float score, Soln prev) {
			for (int k = size() - 1; k > i; k--) {
				solns[k].copy(solns[k - 1]);
			}
			solns[i].setPrevSoln(prev, score);
		}

		protected void add(Entry e, float thisScore) {
			assert (valid);
			if (e == null) {
				add(thisScore);
				return;
			}
			int insertPos = 0;

			for (int i = 0; (i < e.size()) && (insertPos < size()); i++) {

				float score = e.get(i).score + thisScore;
				insertPos = findInsert(insertPos, score, e.get(i));
			}
			// print();
		}

		protected int findInsert(int insertPos, float score, Soln prev) {

			for (; insertPos < size(); insertPos++) {
				System.out.print(score + " " + get(insertPos).score + " | ");
				if (score >= get(insertPos).score) {
					insert(insertPos, score, prev);
					insertPos++;
					break;
				}
			}
			return insertPos;
		}

		protected void add(float thisScore) {
			findInsert(0, thisScore, null);
		}

		protected int numSolns() {
			for (int i = 0; i < solns.length; i++)
				if (solns[i].isClear())
					return i;
			return size();
		}

		public void setValid() {
			valid = true;
		}

		void print() {
			String str = "";
			for (int i = 0; i < size(); i++)
				str += ("[" + i + " " + solns[i].score + " i:" + solns[i].pos
						+ " y:" + solns[i].label + "]");
			System.out.println(str);
		}
	};

	Entry winningLabel[][];
	protected Entry finalSoln;
	protected DoubleMatrix2D Mi;
	protected DoubleMatrix1D Ri;

	void allocateScratch(int numY) {
		Mi = new DenseDoubleMatrix2D(numY, numY);
		Ri = new DenseDoubleMatrix1D(numY);
		winningLabel = new Entry[numY][];
		finalSoln = new Entry(beamsize, 0, 0);
	}

	double fillArray(DataSequence dataSeq, double lambda[], boolean calcScore) {
		double corrScore = 0;
		int numY = model.numY;

		// Mi.assign(0);
		for (int i = 0; i < dataSeq.length(); i++) {
			// compute Mi.

			Trainer.computeLogMi(model.featureGenerator, lambda, dataSeq, i,
					Mi, Ri, false);

			for (int yi = 0; yi < numY; yi++) {
				winningLabel[yi][i].clear();
				winningLabel[yi][i].valid = true;
			}
			for (int yi = model.edgeGen.firstY(i); yi < numY; yi = model.edgeGen
					.nextY(yi, i)) {
				if (i > 0) {
					for (int yp = model.edgeGen.first(yi); yp < numY; yp = model.edgeGen
							.next(yi, yp)) {
						double val = Mi.get(yp, yi) + Ri.get(yi);

						// System.out.println("val "+val);
						winningLabel[yi][i].add(winningLabel[yp][i - 1],
								(float) val);

					}
				} else {
					winningLabel[yi][i].add((float) Ri.get(yi));
				}
			}
			if (calcScore)
				corrScore += (Ri.get(dataSeq.y(i)) + ((i > 0) ? Mi.get(
						dataSeq.y(i - 1), dataSeq.y(i)) : 0));
		}
		return corrScore;
	}

	double fillArray(DataSequence dataSeq, double lambda[], int pos, int label) {
		double corrScore = 0;
		int numY = model.numY;

		// Mi.assign(0);
		for (int i = 0; i < dataSeq.length(); i++) {
			// compute Mi.

			Trainer.computeLogMi(model.featureGenerator, lambda, dataSeq, i,
					Mi, Ri, false);

			for (int yi = 0; yi < numY; yi++) {
				winningLabel[yi][i].clear();
				winningLabel[yi][i].valid = true;
			}
			for (int yi = model.edgeGen.firstY(i); yi < numY; yi = model.edgeGen
					.nextY(yi, i)) {
				if (i > 0) {
					for (int yp = model.edgeGen.first(yi); yp < numY; yp = model.edgeGen
							.next(yi, yp)) {
						double val = Mi.get(yp, yi) + Ri.get(yi);

						// System.out.println("val "+val);
						winningLabel[yi][i].add(winningLabel[yp][i - 1],
								(float) val);
						if ((i == pos) && (yi != label))
							winningLabel[yi][i].solns[0].score = 0;
					}
				} else {
					winningLabel[yi][i].add((float) Ri.get(yi));
				}
			}
			if (false)
				corrScore += (Ri.get(dataSeq.y(i)) + ((i > 0) ? Mi.get(
						dataSeq.y(i - 1), dataSeq.y(i)) : 0));
		}

		return corrScore;
	}
	
	double fillArray(DataSequence dataSeq, double lambda[], int pos, int label, int labelLeft, int labelRight) {
		double corrScore = 0;
		int numY = model.numY;

		// Mi.assign(0);
		for (int i = 0; i < dataSeq.length(); i++) {
			// compute Mi.

			Trainer.computeLogMi(model.featureGenerator, lambda, dataSeq, i,
					Mi, Ri, false);

			for (int yi = 0; yi < numY; yi++) {
				winningLabel[yi][i].clear();
				winningLabel[yi][i].valid = true;
			}
			for (int yi = model.edgeGen.firstY(i); yi < numY; yi = model.edgeGen
					.nextY(yi, i)) {
				if (i > 0) {
					for (int yp = model.edgeGen.first(yi); yp < numY; yp = model.edgeGen
							.next(yi, yp)) {
						double val = Mi.get(yp, yi) + Ri.get(yi);

						// System.out.println("val "+val);
						winningLabel[yi][i].add(winningLabel[yp][i - 1],
								(float) val);
						if ((i == pos) && (yi != label))
							winningLabel[yi][i].solns[0].score = 0;
						if ((i == pos-1) && (yi != labelLeft))
							winningLabel[yi][i].solns[0].score = 0;
						if ((i == pos+1) && (yi != labelRight))
							winningLabel[yi][i].solns[0].score = 0;
					}
				} else {
					winningLabel[yi][i].add((float) Ri.get(yi));
				}
			}
			if (false)
				corrScore += (Ri.get(dataSeq.y(i)) + ((i > 0) ? Mi.get(
						dataSeq.y(i - 1), dataSeq.y(i)) : 0));
		}

		return corrScore;
	}
	
	double fillArray(DataSequence dataSeq, double lambda[], int pos, int pos2, int label, int label2, int flag) {
		double corrScore = 0;
		int numY = model.numY;

		// Mi.assign(0);
		for (int i = 0; i < dataSeq.length(); i++) {
			// compute Mi.

			Trainer.computeLogMi(model.featureGenerator, lambda, dataSeq, i,
					Mi, Ri, false);

			for (int yi = 0; yi < numY; yi++) {
				winningLabel[yi][i].clear();
				winningLabel[yi][i].valid = true;
			}
			for (int yi = model.edgeGen.firstY(i); yi < numY; yi = model.edgeGen
					.nextY(yi, i)) {
				if (i > 0) {
					for (int yp = model.edgeGen.first(yi); yp < numY; yp = model.edgeGen
							.next(yi, yp)) {
						double val = Mi.get(yp, yi) + Ri.get(yi);

						// System.out.println("val "+val);
						winningLabel[yi][i].add(winningLabel[yp][i - 1],
								(float) val);
						if ((i == pos) && (yi != label))
							winningLabel[yi][i].solns[0].score = 0;
						if ((i == pos2) && (yi != label2))
							winningLabel[yi][i].solns[0].score = 0;
					}
				} else {
					winningLabel[yi][i].add((float) Ri.get(yi));
				}
			}
			if (false)
				corrScore += (Ri.get(dataSeq.y(i)) + ((i > 0) ? Mi.get(
						dataSeq.y(i - 1), dataSeq.y(i)) : 0));
		}

		return corrScore;
	}

	protected void setSegment(DataSequence dataSeq, int prevPos, int pos,
			int label) {
		// System.out.println("setting " + pos + " " + label);
		dataSeq.set_y(pos, label);
	}

	public void bestLabelSequence(DataSequence dataSeq, double lambda[]) {
		double corrScore = viterbiSearch(dataSeq, lambda, false);
		assignLabels(dataSeq);
	}

	public void bestLabelSequence(DataSequence dataSeq, double lambda[],
			int pos, int label) {
		double corrScore = viterbiSearch(dataSeq, lambda, pos, label);
		assignLabels(dataSeq);
	}
	
	public void bestLabelSequence(DataSequence dataSeq, double lambda[],
			int pos, int label, int labelLeft, int labelRight) {
		double corrScore = viterbiSearch(dataSeq, lambda, pos, label, labelLeft, labelRight);
		assignLabels(dataSeq);
	}
	
	public void bestLabelSequence(DataSequence dataSeq, double lambda[],
			int pos, int pos2, int label, int label2, int flag) {
		double corrScore = viterbiSearch(dataSeq, lambda, pos, pos2, label, label2, flag);
		assignLabels(dataSeq);
	}

	void assignLabels(DataSequence dataSeq) {
		Soln ybest = finalSoln.get(0);
		ybest = ybest.prevSoln;
		int pos = -1;
		System.out.println("Assigning Labels");
		while (ybest != null) {
			pos = ybest.pos;
			setSegment(dataSeq, ybest.prevPos(), ybest.pos, ybest.label);
			ybest = ybest.prevSoln;
		}
		assert (pos >= 0);
	}

	public double viterbiSearch(DataSequence dataSeq, double lambda[],
			boolean calcCorrectScore) {
		if (Mi == null) {
			allocateScratch(model.numY);
		}
		if ((winningLabel[0] == null)
				|| (winningLabel[0].length < dataSeq.length())) {
			for (int yi = 0; yi < winningLabel.length; yi++) {
				winningLabel[yi] = new Entry[dataSeq.length()];
				for (int l = 0; l < dataSeq.length(); l++)
					winningLabel[yi][l] = new Entry((l == 0) ? 1 : beamsize,
							yi, l);
			}
		}

		double corrScore = fillArray(dataSeq, lambda, calcCorrectScore);
		/*
		 * System.out.println("Winning label\n"); for (int i = 0; i <
		 * winningLabel.length; i++) { for (int j = 0; j < dataSeq.length();
		 * j++) System.out.print(winningLabel[i][j].solns[0].score + " " +
		 * winningLabel[i][j].solns[0].prevLabel() + " " +
		 * winningLabel[i][j].solns[0].label + " | "); System.out.println("\n");
		 * }
		 */
		finalSoln.clear();
		finalSoln.valid = true;
		for (int yi = 0; yi < model.numY; yi++) {
			finalSoln.add(winningLabel[yi][dataSeq.length() - 1], 0);
		}
		/*
		 * System.out.println("final solution"); for (int i = 0; i <
		 * finalSoln.solns.length; i++)
		 * System.out.println(finalSoln.solns[i].label + " " +
		 * finalSoln.solns[i].score + " ");
		 */
		return corrScore;
	}

	public double viterbiSearch(DataSequence dataSeq, double lambda[], int pos,
			int label) {
		if (Mi == null) {
			allocateScratch(model.numY);
		}
		if ((winningLabel[0] == null)
				|| (winningLabel[0].length < dataSeq.length())) {
			for (int yi = 0; yi < winningLabel.length; yi++) {
				winningLabel[yi] = new Entry[dataSeq.length()];
				for (int l = 0; l < dataSeq.length(); l++)
					winningLabel[yi][l] = new Entry((l == 0) ? 1 : beamsize,
							yi, l);
			}
		}

		double corrScore = fillArray(dataSeq, lambda, pos, label);
		System.out.println("Winning label\n");
		for (int i = 0; i < winningLabel.length; i++) {
			for (int j = 0; j < dataSeq.length(); j++)
				System.out.print(winningLabel[i][j].solns[0].score + " "
						+ winningLabel[i][j].solns[0].prevLabel() + " "
						+ winningLabel[i][j].solns[0].label + " | ");
			System.out.println("\n");
		}
		finalSoln.clear();
		finalSoln.valid = true;
		for (int yi = 0; yi < model.numY; yi++) {
			finalSoln.add(winningLabel[yi][dataSeq.length() - 1], 0);
		}
		System.out.println("final solution");
		for (int i = 0; i < finalSoln.solns.length; i++)
			System.out.println(finalSoln.solns[i].label + " "
					+ finalSoln.solns[i].score + " ");
		return corrScore;
	}
	
	public double viterbiSearch(DataSequence dataSeq, double lambda[], int pos,
			int label, int labelLeft, int labelRight) {
		if (Mi == null) {
			allocateScratch(model.numY);
		}
		if ((winningLabel[0] == null)
				|| (winningLabel[0].length < dataSeq.length())) {
			for (int yi = 0; yi < winningLabel.length; yi++) {
				winningLabel[yi] = new Entry[dataSeq.length()];
				for (int l = 0; l < dataSeq.length(); l++)
					winningLabel[yi][l] = new Entry((l == 0) ? 1 : beamsize,
							yi, l);
			}
		}

		double corrScore = fillArray(dataSeq, lambda, pos, label, labelLeft, labelRight);
		System.out.println("Winning label\n");
		for (int i = 0; i < winningLabel.length; i++) {
			for (int j = 0; j < dataSeq.length(); j++)
				System.out.print(winningLabel[i][j].solns[0].score + " "
						+ winningLabel[i][j].solns[0].prevLabel() + " "
						+ winningLabel[i][j].solns[0].label + " | ");
			System.out.println("\n");
		}
		finalSoln.clear();
		finalSoln.valid = true;
		for (int yi = 0; yi < model.numY; yi++) {
			finalSoln.add(winningLabel[yi][dataSeq.length() - 1], 0);
		}
		System.out.println("final solution");
		for (int i = 0; i < finalSoln.solns.length; i++)
			System.out.println(finalSoln.solns[i].label + " "
					+ finalSoln.solns[i].score + " ");
		return corrScore;
	}
	
	public double viterbiSearch(DataSequence dataSeq, double lambda[], int pos, int pos2,
			int label, int label2, int flag) {
		if (Mi == null) {
			allocateScratch(model.numY);
		}
		if ((winningLabel[0] == null)
				|| (winningLabel[0].length < dataSeq.length())) {
			for (int yi = 0; yi < winningLabel.length; yi++) {
				winningLabel[yi] = new Entry[dataSeq.length()];
				for (int l = 0; l < dataSeq.length(); l++)
					winningLabel[yi][l] = new Entry((l == 0) ? 1 : beamsize,
							yi, l);
			}
		}

		double corrScore = fillArray(dataSeq, lambda, pos, pos2, label, label2, flag);
		System.out.println("Winning label\n");
		for (int i = 0; i < winningLabel.length; i++) {
			for (int j = 0; j < dataSeq.length(); j++)
				System.out.print(winningLabel[i][j].solns[0].score + " "
						+ winningLabel[i][j].solns[0].prevLabel() + " "
						+ winningLabel[i][j].solns[0].label + " | ");
			System.out.println("\n");
		}
		finalSoln.clear();
		finalSoln.valid = true;
		for (int yi = 0; yi < model.numY; yi++) {
			finalSoln.add(winningLabel[yi][dataSeq.length() - 1], 0);
		}
		System.out.println("final solution");
		for (int i = 0; i < finalSoln.solns.length; i++)
			System.out.println(finalSoln.solns[i].label + " "
					+ finalSoln.solns[i].score + " ");
		return corrScore;
	}

	int numSolutions() {
		return finalSoln.numSolns();
	}

	Soln getBestSoln(int k) {
		return finalSoln.get(k).prevSoln;
	}
};
