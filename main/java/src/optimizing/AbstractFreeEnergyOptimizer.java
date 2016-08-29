package optimizing;

import java.io.IOException;
import java.util.ArrayList;

import algorithm.MeanFieldForBayesNet;
import bayesNet.BayesNetHandler;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.NumericalDifferentiableFunction;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.SubSequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.data.sequences.SimpleDiscreteSequence;
import de.jstacs.utils.SafeOutputStream;
import models.AbstractPhyloModel;
import models.AlignmentBasedModel;
import projects.dispom.PFMComparator;
import util.EvaluationUtil;

/**
 * This Abstract Optimizer is the basic for all optimizers which should optimize the free energy. It needs a pointer to
 * {@link BayesNetHandler} and to an array of {@link MeanFieldForBayesNet}s. As mentioned it should optimize the free
 * energy, provided by the MeanField implementations.
 * 
 * @author Chaos
 * 
 */
public abstract class AbstractFreeEnergyOptimizer extends NumericalDifferentiableFunction {

    private Sample[] singleSequenceSample;
    public static boolean OPTIMIZE_MEANFIELDS = false;

    /**
     * A pointer to the meanfields which should be optimized. Each element represents one set of positions (e.g. a
     * motif)
     */
    protected MeanFieldForBayesNet[] vlfbn;

    /** A pointer to the inner bayes net of the model */
    protected BayesNetHandler bnh;

    /** A pointer to the model that should be optimized */
    AbstractPhyloModel model;

    /** weights, used to emphasize some positions. If null. All positions are equally weighted */
    protected double[] weights;

    /** the {@link TerminationCondition} is needed to decide when optimizing should be finished. */
    protected TerminationCondition tc;

    protected int dimension;

    public AbstractFreeEnergyOptimizer(double epsilon) throws IllegalArgumentException {
        super(1e-4);
    }

    /**
     * Sets a pointer to the model that should be optimized. And sets a pointer to the underlying
     * {@link BayesNetHandler}. The handler contains the BayesNet, Probability Functions ...
     */
    public void setModel(AbstractPhyloModel model) {
        this.model = model;
        this.bnh = model.getBNH();
    }

    /**
     * Sets a pointer to the weights, which should be used for optimizing
     * 
     * @param val
     */
    public void setWeights(double[] val) {
        this.weights = val;
    }

    /**
     * Sets a pointer to the MeanFields to optimize on.
     * 
     * @param vlfbn
     */
    public void setVariationalLikelihood(MeanFieldForBayesNet[] vlfbn) {
        this.vlfbn = vlfbn;
    }

    /**
     * Calculates the free energy for all MeanFields set and return the sum. The goal should be to minimize this sum.
     * 
     * @return
     * @throws IOException
     */
    protected double getWeightedLogLikelihoodSum() throws EvaluationException {
        double sum = 0;
        for (int u = 0; u < vlfbn.length; u++) {
            try {
                vlfbn[u].initObservation();
                if (OPTIMIZE_MEANFIELDS) {
                    vlfbn[u].optimizeByNormalisation();
                }
            } catch (Exception e) {
                throw new EvaluationException("Observation could not be initialised", e);
            }

            if (weights != null) {
                sum -= vlfbn[u].calcFreeEnergy() * weights[u];
            } else {
                sum -= vlfbn[u].calcFreeEnergy();
            }
        }
        return sum;
    }

    protected SafeOutputStream getOutputStream(boolean output) {
        SafeOutputStream sus;
        if (output) {
            sus = SafeOutputStream.getSafeOutputStream(new SafeOutputStream(System.out));
        } else {
            sus = SafeOutputStream.getSafeOutputStream(new SafeOutputStream(null));
        }
        return sus;
    }

    protected double[][][] getObservableMotifs() throws EvaluationException {
        int species = bnh.getVirtualTree(0).numberOfLeafs;
        // determine PWM per species
        double[][][] observablePWMs = new double[species][][];
        AlignmentBasedModel abm = new AlignmentBasedModel(vlfbn[0].seq.getLength(), (byte) 0,
                vlfbn[0].seq.getAlphabetContainer());
        for (int o = 0; o < species; o++) {
            try {
                abm.train(singleSequenceSample[o], weights);
                observablePWMs[o] = abm.getCondProbs();
            } catch (Exception e) {
                throw new EvaluationException("Could not learn PWM for species " + o, e);
            }
        }
        return observablePWMs;
    }

    protected double[][][] getLearnedMotifs() {
        int species = bnh.getVirtualTree(0).numberOfLeafs;
        // determine learned PWMs
        double[][][] learnedPWMs = new double[species][][];
        for (int o = 0; o < species; o++) {
            learnedPWMs[o] = model.getCondProbsPerSpecies(o);
        }
        return learnedPWMs;
    }

    protected double getNegativeSumDifferenceOfIC() throws EvaluationException {
        int species = bnh.getVirtualTree(0).numberOfLeafs;

        // determine learned PWMs
        double[][][] learnedPWMs = getLearnedMotifs();
        double[][][] observablePWMs = getObservableMotifs();

        double sum = 0;
        for (int o = 0; o < species; o++) {
            sum += Math.abs(EvaluationUtil.determineIC(learnedPWMs[o])
                    - EvaluationUtil.determineIC(observablePWMs[o]));
        }
        return -sum;
    }

    protected double getSumEuclideanDistance() throws EvaluationException {
        int species = bnh.getVirtualTree(0).numberOfLeafs;

        // determine learned PWMs
        double[][][] learnedPWMs = getLearnedMotifs();
        double[][][] observablePWMs = getObservableMotifs();

        PFMComparator.NormalizedEuclideanDistance euclid = new PFMComparator.NormalizedEuclideanDistance();
        double sum = 0;
        for (int o = 0; o < species; o++) {
            sum += euclid.getDistance(learnedPWMs[o], observablePWMs[o], 0);
        }
        return -sum;
    }

    protected double getSumKLDivergences() throws EvaluationException {
        int species = bnh.getVirtualTree(0).numberOfLeafs;

        // determine learned PWMs
        double[][][] learnedPWMs = getLearnedMotifs();
        double[][][] observablePWMs = getObservableMotifs();

        PFMComparator.SymmetricKullbackLeiblerDivergence klMeasure = new PFMComparator.SymmetricKullbackLeiblerDivergence(
                0);
        double sum = 0;
        for (int o = 0; o < species; o++) {
            sum += klMeasure.getDistance(learnedPWMs[o], observablePWMs[o], 0);
        }
        return -sum;
    }

    protected void initSpeciesSamples() throws WrongLengthException, WrongAlphabetException, EmptySampleException {
        int species = bnh.getVirtualTree(0).numberOfLeafs;
        singleSequenceSample = new Sample[species];
        ArrayList<Sequence<?>>[] singleSequencesPerSpescies = new ArrayList[species];
        for (int o = 0; o < species; o++) {
            singleSequencesPerSpescies[o] = new ArrayList<Sequence<?>>();
        }

        for (int u = 0; u < vlfbn.length; u++) {
            // no subsequencing is needed
            // Get parent sequence
            MultiDimensionalDiscreteSequence parent;
            if (vlfbn[u].seq instanceof MultiDimensionalDiscreteSequence) {
                parent = (MultiDimensionalDiscreteSequence) vlfbn[u].seq;
            } else {
                parent = (MultiDimensionalDiscreteSequence) ((SubSequence<?>) vlfbn[u].seq).getParent();
            }

            for (int o = 0; o < species; o++) {
                // generate MultiDimensionalDiscreteSequence that contains only one sequence
                MultiDimensionalDiscreteSequence singleSequenceParent = new MultiDimensionalDiscreteSequence(
                        parent.getAnnotation(), (SimpleDiscreteSequence) parent.getSequence(o));
                // extract subsequence if needed
                if (vlfbn[u].seq instanceof MultiDimensionalDiscreteSequence) {
                    singleSequencesPerSpescies[o].add(singleSequenceParent);
                } else {
                    singleSequencesPerSpescies[o].add(singleSequenceParent.getSubSequence(
                            ((SubSequence<?>) vlfbn[u].seq).getIndex(0), vlfbn[u].seq.getLength()));
                }
            }
        }

        for (int o = 0; o < species; o++) {
            singleSequenceSample[o] = new Sample("Single Sequence Sample for species " + o,
                    singleSequencesPerSpescies[o].toArray(new Sequence<?>[0]));
        }
    }

    /**
     * Overwrites the current TerminationCondition
     * 
     * @param smallDifferenceOfFunctionEvaluationsCondition
     */
    public void setTerminationCondition(TerminationCondition tc) {
        this.tc = tc;
    }

    public abstract void startOptimizing(boolean output);

    @Override
    public int getDimensionOfScope() {
        return this.dimension;
    }
}
