package evolution;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.Normalisation;
import util.Util;

/**
 * Beta representation of the F81 Model. That means that branch lengths represent substitution rates and NOT substitution rates probabilities.
 * 
 * @author mnettling
 *
 */
public class FS81beta extends EvolModel {
    /** the dimenstion of this model */
    private int _dimension;
    /** the underlying stationary distributions */
	private double[][] _pi;
    /** the transition matrix. Calculated by {@link #reinit(double...)} */
	private double[][] _transitionMatrix;

    /**
     * Instantiates a new Felsenstein81 model.
     * 
     * @param dim
     *            the dimension of this model
     */
    protected FS81beta(int dim) {
        _dimension = dim;
        init();
    }

    public void init() {
		_pi = new double[_dimension][_alphabetSize];

        // Initialize with uniform distributions
        for (int i = 0; i < _dimension; i++) {
			Arrays.fill(_pi[i], 1);
			Normalisation.sumNormalisation(_pi[i]);
        }
		_transitionMatrix = new double[_dimension * _alphabetSize][_alphabetSize];
		reinit(0);
    }

    @Override
    public double[][] getUnweightedTransitions() {
        int indexTmp;
		double[][] m = new double[_transitionMatrix.length][_transitionMatrix[0].length];
        for (int i = 0; i < _dimension; i++) {
			for (int a1 = 0; a1 < _alphabetSize; a1++) {
				for (int a2 = 0; a2 < _alphabetSize; a2++) {
					indexTmp = a1 + i * _alphabetSize;
					m[indexTmp][a2] = _pi[i][a2];
                }
            }
        }
        return m;
    }

    @Override
    public int getDimension() {
        return _dimension;
    }

    /**
     * @return a pointer to the underlying transition matrix. The transition matrix might change when calling
     *         {@link #reinit(double...)}
     */
    @Override
    public double[][] getCurrentTransitions() {
		return _transitionMatrix;
    }

    /**
     * @param v
     *            The expected number of mutations per site
     */
    @Override
    public void reinit(double... v) {
        int indexTmp;
        double eps = 1e-4;
        for (int i = 0; i < _dimension; i++) {
			double beta = 1. / (1. - _pi[i][0] * _pi[i][0] - _pi[i][1] * _pi[i][1] - _pi[i][2] * _pi[i][2] - _pi[i][3] * _pi[i][3]);
            double tmpAlpha = Math.exp(-v[0] * beta);
            if (tmpAlpha < 0 + eps) {
                tmpAlpha = 0;
            } else if (tmpAlpha > 1 - eps) {
                tmpAlpha = 1;
            }
			for (int a1 = 0; a1 < _alphabetSize; a1++) {
				for (int a2 = 0; a2 < _alphabetSize; a2++) {
					indexTmp = a1 + i * _alphabetSize;
                    if (a1 == a2) {
						_transitionMatrix[indexTmp][a2] = tmpAlpha + _pi[i][a2] * (1. - tmpAlpha);
                    } else {
						_transitionMatrix[indexTmp][a2] = _pi[i][a2] * (1. - tmpAlpha);
                    }
                }
            }
        }
    }

    @Override
    public double[][] getStatDistr() {
		return _pi;
    }

    @Override
    public void setStatDistr(double[]... A) {
		if (A.length != _pi.length && A[0].length != _pi[0].length) {
            return;
        }
		_pi = Util.arraycopy(A);
		reinit(0);
    }

    /** Writes this FS81Q instance to XML */
    @Override
    public StringBuffer toXML() {
        StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, _pi, "statDistr");
        XMLParser.appendObjectWithTags(xml, _dimension, "dimension");
        return xml;
    }

    /**
     * Generates a FS81P-Instance from the given XML.
     * 
     * @throws NonParsableException
     */
    public static FS81beta getObjectFromXML(StringBuffer xml) throws NonParsableException {
        double[][] stat = (double[][]) XMLParser.extractObjectForTags(xml, "statDistr");
        int dim = (Integer) XMLParser.extractObjectForTags(xml, "dimension");

        FS81beta Otmp = new FS81beta(dim);
        Otmp.setStatDistr(stat);
        return Otmp;
    }

}
