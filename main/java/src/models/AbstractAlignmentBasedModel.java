package models;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.utils.Normalisation;
import util.Util;

public abstract class AbstractAlignmentBasedModel extends PhyloPreparedAbstractModel {

    protected double[][] lnCondProb;
    protected double[][] condProb;
    protected int[] pows;

    /**
     * @param alphabets
     * @param order
     * @param length
     */
    public AbstractAlignmentBasedModel(AlphabetContainer alphabets, byte order, int length) {
        super(alphabets, order, length);
    }

    /**
     * @param xml
     * @throws NonParsableException
     */
    public AbstractAlignmentBasedModel(StringBuffer xml) throws NonParsableException {
        super(xml);
    }

    protected final void initPows() {
        pows = new int[order + 1];
        pows[0] = 1;
        for (int k = 1; k < pows.length; k++) {
            pows[k] = pows[k - 1] * 4;
        }
    }

    /**
     * Returns a pointer to the internal sequenceWeights array.
     * 
     * @return the internal sequenceWeights array.
     */
    public double[] getSequenceWeights() {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    /**
     * Copies the given seqWeights to the internal sequenceWeights array
     * 
     * @param seqWeights
     */
    public void setNewSeqWeights(double[] seqWeights) {
        throw new UnsupportedOperationException("Not implemented yet");
    }

    @Override
    public void setLnCondProb(double[] lambda, int position) {
        if (lambda.length != lnCondProb[position].length) {
            throw new IllegalArgumentException("The given lambda array must be of length "
                    + lnCondProb[position].length + ".");
        } else if (position > length) {
            throw new IllegalArgumentException("The given position must be smaller than the length (" + length
                    + "), but is " + position + ".");
        } else {
            this.lnCondProb[position] = Arrays.copyOf(lambda, lambda.length);
            for (int a = 0; a < lambda.length; a++) {
                this.condProb[position][a] = Math.exp(lambda[a]);
            }
        }
    }

    @Override
    public void setCondProb(double[] condProb, int position) {
        if (condProb.length != this.condProb[position].length) {
            throw new IllegalArgumentException("The given lambda array must be of length "
                    + this.condProb[position].length + ".");
		} else if (position > Math.max(length, order)) {
            throw new IllegalArgumentException("The given position must be smaller than the length (" + length
					+ ") or the order (" + order + ") of the target probabilities, but is " + position + ".");
        } else {
            this.condProb[position] = Arrays.copyOf(condProb, condProb.length);
            for (int a = 0; a < condProb.length; a++) {
                this.lnCondProb[position][a] = Math.log(condProb[a]);
            }
        }
    }

    /** @return a pointer to the underlying conditional probability array at position i */
    @Override
    public double[] getLnCondProb(int i) {
        return lnCondProb[i];
    }

    /** @return a pointer to the underlying conditional probabilities */
    @Override
    public double[][] getLnCondProbs() {
        return lnCondProb;
    }

    /** @return a pointer to the underlying conditional probability array at position i */
    @Override
    public double[] getCondProb(int i) {
        return condProb[i];
    }

    public double[][] getCondProbsPerSpecies(int s) {
        return getCondProbs();
    }

    /** @return a pointer to the underlying conditional probability array at position i */
    @Override
    public double[][] getCondProbs() {
        return condProb;
    }

    @Override
    public void fillLnCondProbRandomly() {
        for (int i = 0; i < condProb.length; i++) {
            for (int k = 0; k < condProb[i].length; k += 4) {
                double[] tmp = Util.getRandomStochVector(4);
                Normalisation.sumNormalisation(tmp);
                for (int l = k; l < k + 4; l++) {
                    condProb[i][l] = tmp[l - k];
                }
            }
            for (int j = 0; j < condProb[i].length; j++) {
                lnCondProb[i][j] = Math.log(condProb[i][j]);
            }
        }
    }
}
