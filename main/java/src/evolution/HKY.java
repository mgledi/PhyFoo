package evolution;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.Normalisation;
import util.Util;

public class HKY extends EvolModel {
    /** the dimenstion of this model */
    private int dimension;
    /** the underlying stationary distributions */
    private double[][] pi;
    /** the transition matrix. Calculated by {@link #reinit(double...)} */
    private double[][] transitionMatrix;

	private double alpha_beta_ratio;

    /**
     * Instantiates a new Felsenstein81 model.
     * 
     * @param dim
     *            the dimension of this model
     */
	protected HKY(int dim, double alpha_beta_ratio) {
        super();
		if (_alphabetSize != 4) {
			throw new IllegalStateException("Only alphabets of size 4 supported by HKY");
		}

		alpha_beta_ratio = alpha_beta_ratio;
        dimension = dim;
        init();
    }

    public void init() {
		pi = new double[dimension][_alphabetSize];
        for (int i = 0; i < dimension; i++) {
			Arrays.fill(pi[i], 1);
			Normalisation.sumNormalisation(pi[i]);
        }
		transitionMatrix = new double[dimension * _alphabetSize][_alphabetSize];
		reinit(1);
    }

    @Override
    public int getDimension() {
        return dimension;
    }

    @Override
    public double[][] getUnweightedTransitions() {
        int indexTmp;
        double[][] m = new double[transitionMatrix.length][transitionMatrix[0].length];
        for (int i = 0; i < dimension; i++) {
			for (int a1 = 0; a1 < _alphabetSize; a1++) {
				for (int a2 = 0; a2 < _alphabetSize; a2++) {
					indexTmp = a1 + i * _alphabetSize;
                    m[indexTmp][a2] = pi[i][a2];
                }
            }
        }
        return m;
    }

    /**
     * @return a pointer to the underlying transition matrix. The transtion matrix might change when calling
     *         {@link #reinit(double...)}
     */
    @Override
    public double[][] getCurrentTransitions() {
        return transitionMatrix;
    }

    @Override
    public void reinit(double... evolDist) {
        double alpha = evolDist[0];
		double beta = alpha * alpha_beta_ratio;


        for (int i = 0; i < dimension; i++) {
			double aPiA = alpha * pi[i][0], aPiC = alpha * pi[i][1], aPiG = alpha * pi[i][2], aPiT = alpha * pi[i][3];
			double bPiA = beta * pi[i][0], bPiC = beta * pi[i][1], bPiG = beta * pi[i][2], bPiT = beta * pi[i][3];

			transitionMatrix[0 + i * _alphabetSize][0] = 1 - (aPiG + bPiT + bPiC); // A->A
			transitionMatrix[0 + i * _alphabetSize][1] = bPiC; // A->C
			transitionMatrix[0 + i * _alphabetSize][2] = aPiG; // A->G
			transitionMatrix[0 + i * _alphabetSize][3] = bPiT; // A->T

			transitionMatrix[1 + i * _alphabetSize][0] = bPiA; // C->A
			transitionMatrix[1 + i * _alphabetSize][1] = 1 - (aPiT + bPiA + bPiG); // C->C
			transitionMatrix[1 + i * _alphabetSize][2] = bPiG; // C->G
			transitionMatrix[1 + i * _alphabetSize][3] = aPiT; // C->T

			transitionMatrix[2 + i * _alphabetSize][0] = aPiA; // G->A
			transitionMatrix[2 + i * _alphabetSize][1] = bPiC; // G->C
			transitionMatrix[2 + i * _alphabetSize][2] = 1 - (aPiA + bPiT + bPiC); // G->G
			transitionMatrix[2 + i * _alphabetSize][3] = bPiT; // G->T

			transitionMatrix[3 + i * _alphabetSize][0] = bPiA; // T->A
			transitionMatrix[3 + i * _alphabetSize][1] = aPiC; // T->C
			transitionMatrix[3 + i * _alphabetSize][2] = bPiG; // T->G
			transitionMatrix[3 + i * _alphabetSize][3] = 1 - (aPiC + bPiA + bPiG); // T->T
        }
    }

    @Override
    public double[][] getStatDistr() {
		return pi;
    }

    @Override
    public void setStatDistr(double[]... A) {
		if (A.length != pi.length && A[0].length != pi[0].length) {
            return;
        }
		pi = Util.arraycopy(A);
		reinit(1);
    }

    @Override
    public StringBuffer toXML() {
        StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, pi, "statDistr");
        XMLParser.appendObjectWithTags(xml, dimension, "dimension");
		XMLParser.appendObjectWithTags(xml, alpha_beta_ratio, "alpha_beta_ratio");
        return xml;
    }

    /**
     * Generates a FS81Q-Instance from the given XML.
     * 
     * @throws NonParsableException
     */
    public static HKY getObjectFromXML(StringBuffer xml) throws NonParsableException {
        double[][] stat = (double[][]) XMLParser.extractObjectForTags(xml, "statDistr");
        int dim = (Integer) XMLParser.extractObjectForTags(xml, "dimension");
		double alpha_beta_ratio = (Integer) XMLParser.extractObjectForTags(xml, "alpha_beta_ratio");

		HKY Otmp = new HKY(dim, alpha_beta_ratio);
        Otmp.setStatDistr(stat);
        return Otmp;
    }
}
