package models;

import java.util.Arrays;

import org.apache.log4j.Logger;

import algorithm.MeanFieldForBayesNet;
import algorithm.SimpleNodeElimination;
import bayesNet.BayesNet;
import bayesNet.BayesNetHandler;
import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.SubSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.io.XMLParser;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.results.NumericalResultSet;
import evolution.EvolModel;
import io.Alphabet;
import io.PhyloSample;
import io.SequenceSpecificDataSelector;
import optimizing.AbstractFreeEnergyOptimizer;
import optimizing.FE_byParameter_ForBayesNodeJstacs;
import optimizing.MotifParamOptimizer;
import optimizing.TemperatureBasedMotifParamOptimizer;
import optimizing.TemperatureOptimizer;
import optimizing.branchlengths.EdgeOptimizer;
import optimizing.branchlengths.GlobalEdgeOptimizer;
import util.MatrixLinearisation;
import util.PWMUtil;
import util.Util;

/**
 * PhylogeneticMotifModel
 */
public class PhyloBayesModel extends AbstractPhyloModel {
	private static Logger LOGGER = Logger.getLogger(PhyloBayesModel.class);

    /** if true, the motif will be trained by calling #train() */
    public static boolean TRAIN_MOTIF = true;

    /**
     * Shows if the edgle-lenghts should be learned during training.<br>
     * 0 = off<br>
     * 1 = train edges for all positions independent<br>
     * 2 = train edges globally (each position in the model will have the same edgelenghts)<br>
     */
    public static int EDGE_LEARNING = 0;

    /** Should the likelihood be calculated instead of free energy (may be very slow) */
    public static boolean CALC_LOGLIKELIHOOD = false;

    /** During the training phase nothing will be done before this iteration is reached */
    public static int SKIP_TRAINING_STEPS = 1;

    /** Shows if cloning is allowed. If this var is setted to false, the clone-method will return <code>this</code> */
    public static boolean CLONING_ALLOWED = true;

    /** which optimizer should be used */
    public static boolean OPTIMIZE_IN_TOTO = false;

    /** If filtered data is assumed another optimizer must be used */
    public static int ASSUME_FILTERED_DATA = 0;

    public static int INTERNAL_TRAINING_REPEATS = 1;
    public static int INTERNAL_TRAINING_REPEATS_FIRST_STEP = 1;

    private boolean trained = true;

    protected int trainingStep = 0;

    /** Remember the last weights for training */
    private double[] lastWeights;

    /** Remember the last Sample for training */
    private Sample lastData;

    /**
     * Constructor for instantiating a {@link PhyloBayesModel}.
     * 
     * @param motifLength
     *            the length of the motif
     * @param order
     *            the order of the motif
     * @param newickString
     *            the newick representation of the phylogenetic tree to use
     * @param con
     *            the underlying AlphabetContainer
     */
    public PhyloBayesModel(int motifLength, byte order, String newickString, AlphabetContainer con) {
        super(con, order, motifLength);
        bnh = buildStructureSimpleForBasicModel(newickString);
        fillParameter(true);
    }

    /** Constructor to instantiate a {@link PhyloBayesModel} from its String representation */
    public PhyloBayesModel(StringBuffer xml) throws Exception {
        super(xml);
    }

    /** @return PhyloBayesModel#getLength() */
    public int getMotifLength() {
        return getLength();
    }

    /**
     * This method trains the Model on the last known data and the last known weights. That means, that you clearly have
     * to run at least on EM-Step with {@link SingleHiddenMotifMixture}.
     * 
     * @throws CloneNotSupportedException
     */
    public void trainOnLastData() throws CloneNotSupportedException {
        this.train(lastData, lastWeights);
    }

    // ######################## Methoden um das Model in Jstacs inzubinden
    // ################################ */

    public void train(Sample data, double[] weights) throws CloneNotSupportedException {
        // ##################### initialize Training ##################################################
        long startTime = System.currentTimeMillis();
        lastWeights = Arrays.copyOf(weights, weights.length);
        lastData = data;
        if (trainingStep < SKIP_TRAINING_STEPS) {
			LOGGER.info("Training ignored in training step " + trainingStep);
            trainingStep++;
            return;
        }

        SequenceSpecificDataSelector tdm = new SequenceSpecificDataSelector(data, weights);
        tdm.prepareForTraining(PROB_THRESH_FOR_WEIGHTS, MOTIF_QUALITY_THRESHOLD);
		LOGGER.info("Using " + (tdm.weightsForTraining.length) + " of " + weights.length + " for training.");
        MeanFieldForBayesNet[] tmpPointer = new MeanFieldForBayesNet[tdm.weightsForTraining.length];
        double[] myWeights = tdm.weightsForTraining;

        // ##################### start Training ##################################################
        // optimize markov model parameters
        for (int s = 0; s < INTERNAL_TRAINING_REPEATS || trainingStep == 0 && s < INTERNAL_TRAINING_REPEATS_FIRST_STEP; s++) {
            for (int m = 0; m < tdm.weightsForTraining.length; m++) {
                tmpPointer[m] = this.getMFfromCache(tdm.seqsForTraining[m]);
                tmpPointer[m].initObservation();
                tmpPointer[m].optimizeByNormalisation();
                myWeights[m] = tdm.weightsForTraining[m];
            }
			LOGGER.info("Optimizing " + tmpPointer.length + " MeanFieldForBayesNet for PhyloBayesModel.");

            if (TRAIN_MOTIF) {
                if (OPTIMIZE_IN_TOTO) {
                    AbstractFreeEnergyOptimizer optPar;
                    if (ASSUME_FILTERED_DATA > 0) {
                        if (ASSUME_FILTERED_DATA != 1) {
                            throw new IllegalArgumentException("Cant learn local temperatures and motif at once.");
                        }
						LOGGER.info("Use TemperatureBasedMotifParamOptimizer.");
                        optPar = new TemperatureBasedMotifParamOptimizer();
                    } else {
						LOGGER.info("Use MotifParamOptimizer.");
                        optPar = new MotifParamOptimizer();
                    }
                    optPar.setVariationalLikelihood(tmpPointer);
                    optPar.setModel(this);
                    optPar.setWeights(myWeights);
                    optPar.startOptimizing(false);
                } else {
                    AbstractFreeEnergyOptimizer optPar;
					LOGGER.info("Use FE_byParameter_ForBayesNodeJstacs.");
                    optPar = new FE_byParameter_ForBayesNodeJstacs();
                    optPar.setVariationalLikelihood(tmpPointer);
                    optPar.setModel(this);
                    optPar.setWeights(myWeights);

                    int order[] = PWMUtil.getRandomOrder(length);
                    for (int k = 0; k < bnh.motifLength; k++) {
                        ((FE_byParameter_ForBayesNodeJstacs) optPar).initPosition(order[k]);
                        optPar.startOptimizing(false);
                    }

                    if (ASSUME_FILTERED_DATA > 0) {
						LOGGER.info("Use TemperatureOptimizer.");
                        optPar = new TemperatureOptimizer(ASSUME_FILTERED_DATA == 2);
                        optPar.setVariationalLikelihood(tmpPointer);
                        optPar.setModel(this);
                        optPar.setWeights(myWeights);
                        optPar.startOptimizing(false);
                    }
                }
            }

            // optimize edges
            if (EDGE_LEARNING == 1) {
                EdgeOptimizer edgeOptimizer = new EdgeOptimizer(false, EDGE_LEARNING_ONE_LENGTH);
                edgeOptimizer.setModel(this);
                edgeOptimizer.setVariationalLikelihood(tmpPointer);
                edgeOptimizer.setWeights(myWeights);
                for (int k = 0; k < bnh.motifLength; k++) {
                    edgeOptimizer.initPosition(k);
                    edgeOptimizer.startOptimizing(false);
                }
            } else if (EDGE_LEARNING == 2) {
                // TODO: GlobalSingleEdgeOptimizer
                GlobalEdgeOptimizer edgeOptimizer = new GlobalEdgeOptimizer();
                edgeOptimizer.setModel(this);
                edgeOptimizer.setVariationalLikelihood(tmpPointer);
                edgeOptimizer.setWeights(myWeights);
                edgeOptimizer.startOptimizing(false);
            }
        }
        if (EDGE_LEARNING > 0) {
			LOGGER.info("----- New NewickStrings -----");
            for (int k = 0; k < bnh.motifLength; k++) {
				LOGGER.info(bnh.getVirtualTree(k).getNewickString());
            }
        }
		LOGGER.info("Finished Training in " + (System.currentTimeMillis() - startTime) + " milli seconds.");
        trainingStep++;
    }

    @Override
    public double getLogProbFor(Sequence sequence, int startpos, int endpos)
            throws IllegalArgumentException, NotTrainedException {
        MeanFieldForBayesNet tmpPointer;
        this.check(sequence, startpos, endpos);
        Sequence<int[]> subSequence = sequence.getSubSequence(startpos, endpos - startpos + 1);
        MultiDimensionalDiscreteSequence parent;
        if (sequence instanceof MultiDimensionalDiscreteSequence) {
            parent = (MultiDimensionalDiscreteSequence) sequence;
        } else {
            parent = (MultiDimensionalDiscreteSequence) ((SubSequence<?>) subSequence).getParent();
        }
        String[] oldTopologies = null;
        String[] newTopologies = null;
        // if a topology from the dataset is available, use this to calculate the likelihood
        if (parent.learnedTopology != null) {
            oldTopologies = new String[this.getLength()];
            newTopologies = new String[this.getLength()];

            for (int i = 0; i < oldTopologies.length; i++) {
                oldTopologies[i] = getBNH().getVirtualTree(i).getNewickString();
            }
            for (int i = 0; i < oldTopologies.length; i++) {
                newTopologies[i] = PhyloSample.combineBranchLengths(parent.learnedTopology, oldTopologies[i], 1.0);
            }
            this.reinitEdgelengths(newTopologies);
            // this.reinitEdgelengths(parent.learnedTopology );
        }

        tmpPointer = this.getMFfromCache(subSequence);
        tmpPointer.initObservation();

        double logLikelihood = 0;
        if (CALC_LOGLIKELIHOOD) {
            SimpleNodeElimination.instantiate(bnh.getNet());
            logLikelihood = tmpPointer.calcLogLikelihood();
        } else {
            tmpPointer.optimizeByNormalisation();
            logLikelihood = -tmpPointer.calcFreeEnergy();
        }

        // reset toplogies, if topologies from data set were uesed
        if (parent.learnedTopology != null) {
            this.reinitEdgelengths(oldTopologies);
        }
        return logLikelihood;
    }

    @Override
    public double getLogPriorTerm() throws Exception {
        return 0;
    }

    /** returns the name of the Instance */
    public String getInstanceName() {
        return "PhyloBayesModel";
    }

    /** true, if the Model is trained */
    @Override
    public boolean isTrained() {
        return trained;
    }

    @Override
    public NumericalResultSet getNumericalCharacteristics() throws Exception {
        throw new UnsupportedOperationException("This method is not implemented yet.");
    }

    @Override
    public StringBuffer toXML() {
        StringBuffer xml = bnh.toXML();
        XMLParser.addTags(xml, "_bnh");
        XMLParser.appendObjectWithTags(xml, order, "_order");
        XMLParser.appendObjectWithTags(xml, length, "modelLength");
        XMLParser.appendObjectWithTags(xml, alphabets, "alphabets");
        return xml;
    }

    @Override
    protected void fromXML(StringBuffer xml) throws NonParsableException {
        length = (Integer) XMLParser.extractObjectForTags(xml, "modelLength");
        alphabets = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "alphabets");
        order = (Byte) XMLParser.extractObjectForTags(xml, "_order");

        bnh = new BayesNetHandler(XMLParser.extractForTag(xml, "_bnh"));
        bnh.setSimpleRoot();
        trained = true;
    }

    @Override
    public PhyloBayesModel clone() {
        if (!CLONING_ALLOWED)
            return this;
        StringBuffer xml = this.toXML();
        try {
            PhyloBayesModel pbm = new PhyloBayesModel(xml);
            return pbm;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return this;
    }

    // #########################################################################
    // ########################### Hilfsfunktionen #############################
    // #########################################################################

    /** fills the underlying bayes net with random or uniformly distributed parameters */
    public void fillParameter(boolean randomly) {
        if (randomly) {
			LOGGER.info("Setting parameters randomly for PhyloBayesModel.");
        } else {
			LOGGER.info("Setting parameters equal distributed");
        }
        EvolModel gfs = EvolModel.getInstance(1);
        for (int i = 0; i < bnh.getVirtualTrees().length; i++) {
            gfs = EvolModel
                    .getInstance((int) Math.pow(Alphabet.size, bnh.getVirtualTree(i).getNode(0).numberOfParents));
            bnh.getVirtualTree(i).setEvolModel(gfs);

            if (randomly) { // set random parameters
                bnh.getVirtualTree(i).getEvolModel()
                        .setStatDistr(Util.getRandomStochMatrix(gfs.getDimension(), Alphabet.size));
            } else { // set equal distributed parameters
                bnh.getVirtualTree(i).getEvolModel()
                        .setStatDistr(Util.getEqualStochMatrix(gfs.getDimension(), Alphabet.size));
            }
            bnh.getVirtualTree(i).initParametersFromGF();
        }
    }

    /**
     * build graphical structure, using at each position the given newickstring
     */
    private BayesNetHandler buildStructureSimpleForBasicModel(String newick) {
        BayesNet wholeNet = new BayesNet("PhyloBayes");
        BayesNetHandler bnh = new BayesNetHandler(wholeNet);

        for (int i = 0; i < length; i++) {
            bnh.addBayesNet(io.NewickToBayesNet.getTree(newick, "pos_" + i), newick);
            // connect with predecessor
            if (i > 0 && order == 1) {
                bnh.connectVirtualTrees(bnh.getVirtualTree(i - 1), bnh.getVirtualTree(i));
            }
        }
        // bnh.topologicalSort();
        bnh.setSimpleRoot();
        return bnh;
    }

    /** Draw random motif, with the underlying parameters */
    public String[] emitSample() {
        return emitSample(null);
    }

    /** Draw random motif, with the underlying parameters using the given old observation for the underlying bayesnet */
    public String[] emitSample(int[] observation) {
        int org = bnh.getVirtualTree(0).numberOfLeafs;
        String[] seq = new String[org];
        for (int o = 0; o < org; o++) {
            seq[o] = "";
        }
        // generate full observation
        int[] fullObs = null;
        if (observation == null) {
            fullObs = bnh.drawFullObservation();
        } else {
            fullObs = bnh.drawFullObservation(observation);
        }
        for (int i = 0; i < length; i++) {
            for (int o = 0; o < org; o++) {
				seq[o] += alphabets.getSymbol(0, fullObs[bnh.getVirtualTree(i).getLeaf(o).nodeNumber]);
            }
        }
        return seq;
    }

    public void setTemperature(double t) {
        for (int i = 0; i < length; i++) {
            bnh.getVirtualTree(i).TEMPERATURE = t;
        }
    }

    /* (non-Javadoc)
     * 
     * @see de.jstacs.models.discrete.DiscreteGraphicalModel#check(de.jstacs.data .Sequence, int, int) */
    protected void check(Sequence<?> sequence, int startpos, int endpos) throws NotTrainedException,
            IllegalArgumentException {
        if (!trained) {
            throw new NotTrainedException();
        } else if (!alphabets.checkConsistency(sequence.getAlphabetContainer().getSubContainer(startpos,
                endpos - startpos + 1))) {
            throw new IllegalArgumentException("This sequence is not possible with the given alphabet.");
        } else if (startpos < 0) {
            throw new IllegalArgumentException("This startposition is impossible. Try: 0 <= startposition");
        } else if (startpos > endpos || endpos >= sequence.getLength()) {
            throw new IllegalArgumentException(
                    "This endposition is impossible. Try: startposition <= endposition < sequence.length");
        } else if (endpos - startpos + 1 != length) {
            throw new IllegalArgumentException("Expected length of sequence is " + length + ". But was "
                    + (endpos - startpos + 1));
        }
    }

	/** Returns a array of probability-arrays */
    public double[][] getProbabilites() {
        double[][] m = new double[this.length][];
        for (int i = 0; i < m.length; i++) {
            double[] tmp = MatrixLinearisation.linearize(bnh.getVirtualTree(i).getEvolModel().getStatDistr());
            m[i] = Arrays.copyOf(tmp, tmp.length);
        }
        return m;
    }

    public double[] getTemperatures() {
        double[] temps = new double[length];
        for (int i = 0; i < temps.length; i++) {
            temps[i] = Util.round(getBNH().getVirtualTree(i).TEMPERATURE, 3);
        }
        return temps;
    }
}
