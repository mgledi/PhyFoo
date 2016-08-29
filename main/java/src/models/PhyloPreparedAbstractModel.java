package models;

import java.util.Arrays;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.models.AbstractModel;
import de.jstacs.models.discrete.homogeneous.HomogeneousMM;
import de.jstacs.models.discrete.homogeneous.parameters.HomMMParameterSet;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.models.discrete.inhomogeneous.parameters.BayesianNetworkModelParameterSet;

/**
 */
public abstract class PhyloPreparedAbstractModel extends AbstractModel {

	public static double MOTIF_QUALITY_THRESHOLD = 0.0;
	public static double PROB_THRESH_FOR_WEIGHTS = 0.8;

    /** The maximal markov order. */
    protected byte order;

    public PhyloPreparedAbstractModel(AlphabetContainer alphabets, byte order, int length) {
        super(alphabets, length);
        this.order = order;
    }

    public PhyloPreparedAbstractModel(StringBuffer xml) throws NonParsableException {
        super(xml);
    }

    /**
     * This method converts this {@link PhyloPreparedAbstractModel} to a {@link BayesianNetworkModel} or a
     * {@link HomogeneousMM}.
     * 
     * @param initTrained
     * @param expectedAlphabet
     * @return
     * @throws Exception
     */
    public AbstractModel toJstacsModel(Sample initTrained, AlphabetContainer expectedAlphabet) throws Exception {
        double[][] motif = this.getCondProbs();
        double[][] lnmotif = this.getLnCondProbs();

        if (this.getLength() > 0) {
            BayesianNetworkModelParameterSet params = new BayesianNetworkModelParameterSet(
                    expectedAlphabet,
                    this.getLength(),
                    1,
                    "foreground model",
                    ModelType.IMM, (byte)
                    this.order,
                    LearningType.ML_OR_MAP);
            BayesianNetworkModel model = new BayesianNetworkModel(params);
            model.train(initTrained.getInfixSample(0, model.getLength())); // must be done to initialize constraints
            for (int i = 0; i < this.getLength(); i++) {
                model.constraints[i].freq = Arrays.copyOf(motif[i], motif[i].length);
                model.constraints[i].lnFreq = Arrays.copyOf(lnmotif[i], lnmotif[i].length);
            }
            return model;
        } else {
            HomogeneousMM model = new HomogeneousMM(new HomMMParameterSet(expectedAlphabet, 1, null, this.getOrder()));
            model.train(initTrained);
            for (int i = 0; i <= this.getOrder(); i++) {
                model.condProb[i].freq = Arrays.copyOf(motif[i], motif[i].length);
                model.condProb[i].lnFreq = Arrays.copyOf(lnmotif[i], lnmotif[i].length);
            }
            return model;
        }
    }

    /**
     * Sets the logarithmic conditional probability for position i
     * 
     * @param lambda
     *            the lambda vector to set. The length of the lambda vector is always a multiplier of the
     *            alphabet-length
     * @param position
     *            the position of the given lambda vector
     */
    public abstract void setLnCondProb(double[] lambda, int position);

    /**
     * Sets the logarithmic conditional probabilities
     * 
     * @param lambdas
     *            the lambda vector to set. The length of the lambda vector is always a multiplier of the
     *            alphabet-length
     */
    public void setLnCondProbs(double[][] lambdas) {
        for (int i = 0; i < lambdas.length; i++) {
            this.setLnCondProb(lambdas[i], i);
        }
    }

    @Override
    public byte getMaximalMarkovOrder() {
        return getOrder();
    }

    /**
     * @return the order of the model
     */
    public byte getOrder() {
        return order;
    }

    /** @return a the underlying conditional probability array at position i */
    public abstract double[][] getLnCondProbs();

    /** @return a the underlying conditional logarithmic probability array at position i */
    public abstract double[] getLnCondProb(int i);

    /**
     * Sets the conditional probability for position i
     * 
     * @param condProb
     *            the lambda vector to set. The length of the lambda vector is always a multiplier of the
     *            alphabet-length
     * @param position
     *            the position of the given lambda vector
     */
    public abstract void setCondProb(double[] condProb, int position);

    /**
     * Sets the conditional probability for position i
     * 
     * @param condProbs
     *            the lambda vector to set. The length of the lambda vector is always a multiplier of the
     *            alphabet-length
     */
    public void setCondProbs(double[][] condProbs) {
        for (int i = 0; i < condProbs.length; i++) {
            this.setCondProb(condProbs[i], i);
        }
    }

    /** @return a copy of the underlying conditional probability array at position i */
    public abstract double[] getCondProb(int i);

    /** @return a copy of the underlying conditional probability array at position i */
    public abstract double[][] getCondProbs();

    public abstract double[][] getCondProbsPerSpecies(int s);

    public abstract void fillLnCondProbRandomly();
}
