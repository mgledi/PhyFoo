package evolution;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import de.jstacs.utils.Normalisation;
import util.Util;

/**
 * Alpha representation of the F81 Model. That means that branch lengths represent substitution probabilities and NOT substitution rates.
 * 
 * @author Martin Nettling
 *
 */
public class FS81alpha extends EvolModel {
	/** the dimension of this model */
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
	protected FS81alpha(int dim) {
        super();
		_dimension = dim;
        init();
    }

    public void init() {
		_pi = new double[_dimension][_alphabetSize];
		for (int i = 0; i < _dimension; i++) {
			Arrays.fill(_pi[i], 1);
			Normalisation.sumNormalisation(_pi[i]);
        }
		_transitionMatrix = new double[_dimension * _alphabetSize][_alphabetSize];
		reinit(1);
    }

    @Override
    public int getDimension() {
		return _dimension;
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

    /**
     * @return a pointer to the underlying transition matrix. The transtion matrix might change when calling
     *         {@link #reinit(double...)}
     */
    @Override
    public double[][] getCurrentTransitions() {
		return _transitionMatrix;
    }

    @Override
    public void reinit(double... evolDist) {
        int indexTmp;
        double alpha = evolDist[0];

		for (int i = 0; i < _dimension; i++) {
			for (int a1 = 0; a1 < _alphabetSize; a1++) {
				for (int a2 = 0; a2 < _alphabetSize; a2++) {
					indexTmp = a1 + i * _alphabetSize;
                    if (a1 == a2) {
						_transitionMatrix[indexTmp][a2] = 1 - alpha + _pi[i][a2] * alpha;
                    } else {
						_transitionMatrix[indexTmp][a2] = _pi[i][a2] * alpha;
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
		reinit(1);
    }

    @Override
    public StringBuffer toXML() {
        StringBuffer xml = new StringBuffer();
		XMLParser.appendObjectWithTags(xml, _pi, "statDistr");
		XMLParser.appendObjectWithTags(xml, _dimension, "dimension");
        return xml;
    }

    /**
	 * Generates a F81Q-Instance from the given XML.
	 * 
	 * @throws NonParsableException
	 */
    public static FS81alpha getObjectFromXML(StringBuffer xml) throws NonParsableException {
        double[][] stat = (double[][]) XMLParser.extractObjectForTags(xml, "statDistr");
        int dim = (Integer) XMLParser.extractObjectForTags(xml, "dimension");

		FS81alpha Otmp = new FS81alpha(dim);
        Otmp.setStatDistr(stat);
        return Otmp;
    }
}
