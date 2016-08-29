package models;

import java.util.Arrays;

import org.apache.log4j.Logger;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.RecursiveSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.Normalisation;
import io.SequenceSpecificDataSelector;

/**
 * This Model is the optimized version of {@link PhyloBayesModel} with branch lengths equal to 1.
 * 
 * @author Martin Nettling
 * 
 */
public class AlignmentBasedModel extends AbstractAlignmentBasedModel {
	private static Logger LOGGER = Logger.getLogger(PhyloBayesModel.class);

    public static boolean IGNORE_GAPS_IN_TRAINING = false;
    public static boolean IGNORE_GAPS_IN_LIKELIHOOD_CALC = false;

    public static int SKIP_TRAINING_STEPS = 1;
	public static double PSEUDO_COUNT = 1e-8;

    private int[] ACCESS_TABLE;
    int trainingStep;
	// only used for hack things
	public double[] lastWeights;
	public Sample lastSequences;

    /**
     * Instanitates a new {@link AlignmentBasedModel} and initialises its underlying conditional probabilities randomly.
     * 
     * @param length
     * @param order
     * @param alphabet
     */
    public AlignmentBasedModel(int length, byte order, AlphabetContainer alphabet) {
        super(alphabet, order, length);
        this.trainingStep = 0;
        this.order = order;
        initPows();

        lnCondProb = new double[length][];
        condProb = new double[length][];

        int m = 0;
        for (; m < order; m++) {
            lnCondProb[m] = new double[pows[m] * 4];
            condProb[m] = new double[pows[m] * 4];
        }
        for (; m < length; m++) {
            lnCondProb[m] = new double[pows[order] * 4];
            condProb[m] = new double[pows[order] * 4];
        }
        fillLnCondProbRandomly();
        int size = 0;
        for (int i = 0; i < pows.length; i++) {
            size += pows[i] * 4;
        }
        ACCESS_TABLE = (new int[size]);
    }

    public AlignmentBasedModel(StringBuffer xml) throws NonParsableException {
        super(xml);
        int size = 0;
        for (int i = 0; i < pows.length; i++) {
            size += pows[i] * 4;
        }
        ACCESS_TABLE = (new int[size]);
    }

    @Override
    public void train(Sample data, double[] weights) throws Exception {
		this.lastWeights = Arrays.copyOf(weights, weights.length);
		this.lastSequences = data;
        if (trainingStep++ < SKIP_TRAINING_STEPS) {
			LOGGER.info("Skipping training step " + trainingStep);
            return;
        }
        // prepare data
        SequenceSpecificDataSelector tdm = new SequenceSpecificDataSelector(data, weights);
        tdm.prepareForTraining(PROB_THRESH_FOR_WEIGHTS, MOTIF_QUALITY_THRESHOLD);
        // System.out.println("Using " + (tdm.weightsForTraining.length) + " of " + data.getNumberOfElements()
        // + " for training.");

        double[/* motif-position */][/* base pair-hash */] counts = new double[length][];
        for (int i = 0; i < length; i++) {
            counts[i] = new double[lnCondProb[i].length];
            Arrays.fill(counts[i], PSEUDO_COUNT);
        }

        Sequence<?>[] sequences = tdm.seqsForTraining;
        double[] myWeights = tdm.weightsForTraining;
		LOGGER.info("Using " + (tdm.weightsForTraining.length) + " of " + weights.length + " for training.");

        for (int w = 0; w < sequences.length; w++) {
            MultiDimensionalDiscreteSequence parent;
            int start = 0, end = length;

            if (sequences[w] instanceof MultiDimensionalDiscreteSequence) {
                parent = (MultiDimensionalDiscreteSequence) sequences[w];
                end = parent.getLength() - 1;
            } else {
                RecursiveSequence<?> subseq = (RecursiveSequence<?>) sequences[w];
                parent = (MultiDimensionalDiscreteSequence) subseq.getParent();
                start = subseq.getIndex(0);
                end = subseq.getIndex(length - 1);
            }

            // run over all positions in the given sequence (length)
            for (int o = 0; o < parent.getNumberOfSequences(); o++) {
                for (int l = start; l <= end; l++) {
                    // k,o is the position in the sequence -> determine count hash
                    Arrays.fill(ACCESS_TABLE, 0);
                    int accessHead = 0, preHead = 0;
                    // observations.add(new int[order + 1]);
                    // Arrays.fill(observations.get(0), -1); // at least one observation can be made
                    for (int k = 0; k <= order && k <= (l - start); k++) {
                        if (parent.containerForPhyloBayes[l - k][o] == 4 && IGNORE_GAPS_IN_TRAINING) {
                            // System.out.println("ignoring shit");
                            accessHead = -1; // skip summation if IGNORE_GAPS_IN_TRAINING was true
                            break;
                        } else if (parent.containerForPhyloBayes[l - k][o] == 4) {
                            int oldPreHead = preHead;
                            preHead = accessHead;
                            // here is an error: when adding new observations, old ones must be removed
                            for (int i = oldPreHead; i < preHead || preHead == 0 && i == 0; i++) {
                                for (int a = 0; a < 4; a++) {
                                    ACCESS_TABLE[accessHead++] = ACCESS_TABLE[i] + a * pows[k];
                                }
                            }
                        } else {  // add to all observations the current base
                            accessHead = accessHead == 0 ? 1 : accessHead; // if no observation was made yet, add a observation
                            for (int i = preHead; i < accessHead; i++) {
                                ACCESS_TABLE[i] += parent.containerForPhyloBayes[l - k][o] * pows[k];
                            }
                        }
                    }
                    // increment count-hash by the weight[w]
                    for (int i = preHead; i < accessHead; i++) {
                        counts[l - start][ACCESS_TABLE[i]] += (myWeights == null ? 1. : myWeights[w]) / (accessHead - preHead);
                    }
                }
            }
        }
        // Normalization
        for (int i = 0; i < length; i++) {
            for (int k = 0; k < counts[i].length; k += 4) {
                double[] tmp = Arrays.copyOfRange(counts[i], k, k + 4);
                Normalisation.sumNormalisation(tmp);
                for (int l = k; l < k + 4; l++) {
                    condProb[i][l] = tmp[l - k];
                }
            }
            for (int j = 0; j < counts[i].length; j++) {
                lnCondProb[i][j] = Math.log(condProb[i][j]);
            }
        }
    }

    /**
     * Calculates the log-probability of occurance of the motif at position posSeq in the given sequence.
     */
    private double getLnCondProb(Sequence<?> sequence, int posSeq) {
        double sum = 0;
        for (int i = 0; i < length; i++) {
            sum += getLnCondProb(sequence, posSeq + i, i);
        }
        return sum;
    }

    /**
     * Estimates the conditional probability for nucleotide at position posSeq in the given sequence. The
     * log-probability from posMot in the motif is returned.
     * 
     * @return
     */
    private double getLnCondProb(Sequence<?> sequence, final int posSeq, int posMot) {
        if (AlignmentBasedModel.isCompletylUnObserved(sequence, posSeq - Math.min(order, posMot), posSeq)) {
            return Math.log(0.25);
        }

        // prepare all int[] - arrays, which are represented by the subsequence at posSeq
        Arrays.fill(ACCESS_TABLE, 0);
        int accessHead = 0, preHead = 0; // initially there is no observation
        for (int k = 0; k <= order && k <= posMot; k++) {
            if (sequence.discreteVal(posSeq - k) == 4 && IGNORE_GAPS_IN_LIKELIHOOD_CALC) {
                return 0; // HACKY
            } else if (sequence.discreteVal(posSeq - k) == 4) {
                int oldPreHead = preHead;
                preHead = accessHead;
                // add to all yet made observations all four bases 
                for (int i = oldPreHead; i < preHead || preHead == 0 && i == 0; i++) {
                    for (int a = 0; a < 4; a++) {
                        ACCESS_TABLE[accessHead++] = ACCESS_TABLE[i] + a * pows[k];
                    }
                }
            } else {  // add to all observations the current base
                accessHead = accessHead == 0 ? 1 : accessHead; // if no observation was made yet, add a observation
                for (int i = preHead; i < accessHead ; i++) {
                    ACCESS_TABLE[i] += (sequence.discreteVal(posSeq - k) * pows[k]);
                }
            }
        }

        double sumProb = 0;
        for (int i = preHead; i < accessHead; i++) {
            sumProb += condProb[posMot][ACCESS_TABLE[i]];
        }
        return Math.log(sumProb / (accessHead - preHead));
    }

    /* (non-Javadoc)
     * 
     * @see de.jstacs.models.AbstractModel#getLogProbFor(de.jstacs.data.Sequence, int, int) */
    @Override
    public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws NotTrainedException, Exception {
        check(sequence, startpos, endpos);
        MultiDimensionalDiscreteSequence parent;

        if (sequence instanceof MultiDimensionalDiscreteSequence) {
            parent = (MultiDimensionalDiscreteSequence) sequence;
        } else {
            RecursiveSequence<?> subseq = (RecursiveSequence<?>) sequence;
            parent = (MultiDimensionalDiscreteSequence) subseq.getParent();
            startpos = subseq.getIndex(startpos);
            endpos = subseq.getIndex(endpos);
        }

        double logSum = 0;
//        if(startpos==0) System.out.println("LogCalc " + Arrays.toString(sequenceWeights));
        for (int o = 0; o < parent.getNumberOfSequences(); o++) {
            logSum += this.getLnCondProb(parent.getSequence(o), startpos);
        }

        // sum up all nucleotides and probabilities
        return logSum;
    }

    @Override
    public String getInstanceName() {
        return AlignmentBasedModel.class.getSimpleName();
    }

    @Override
    public StringBuffer toXML() {
        StringBuffer xml = new StringBuffer();
        XMLParser.appendObjectWithTags(xml, order, "order");
        XMLParser.appendObjectWithTags(xml, length, "modelLength");
        XMLParser.appendObjectWithTags(xml, alphabets, "alphabets");
        XMLParser.appendObjectWithTags(xml, lnCondProb, "lnCondProb");
        XMLParser.appendObjectWithTags(xml, condProb, "condProb");
        XMLParser.appendObjectWithTags(xml, pows, "pows");
        return xml;
    }

    @Override
    protected void fromXML(StringBuffer xml) throws NonParsableException {
        order = (Byte) XMLParser.extractObjectForTags(xml, "order");
        length = (Integer) XMLParser.extractObjectForTags(xml, "modelLength");
        alphabets = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "alphabets");
        lnCondProb = (double[][]) XMLParser.extractObjectForTags(xml, "lnCondProb");
        condProb = (double[][]) XMLParser.extractObjectForTags(xml, "condProb");
        pows = (int[]) XMLParser.extractObjectForTags(xml, "pows");
    }

    @Override
    public double getLogPriorTerm() throws Exception {
        return 0;
    }

    @Override
    public boolean isTrained() {
        return true;
    }

    @Override
    public NumericalResultSet getNumericalCharacteristics() throws Exception {
        return null;
    }

    @Override
    public AlignmentBasedModel clone() {
        StringBuffer xml = this.toXML();
        try {
            AlignmentBasedModel abm = new AlignmentBasedModel(xml);
            return abm;
        } catch (Exception e) {
            e.printStackTrace();
        }
        return this;
    }

    /* (non-Javadoc)
     * 
     * @see de.jstacs.models.discrete.DiscreteGraphicalModel#check(de.jstacs.data .Sequence, int, int) */
    protected void check(Sequence<?> sequence, int startpos, int endpos) throws IllegalArgumentException {
        if (!alphabets.checkConsistency(sequence.getAlphabetContainer().getSubContainer(startpos,
                endpos - startpos + 1))) {
            throw new IllegalArgumentException("This sequence is not possible with the given alphabet.");
        } else if (startpos < 0) {
            throw new IllegalArgumentException("This startposition is impossible. Try: 0 <= startposition");
        } else if (startpos > endpos || endpos >= sequence.getLength()) {
            throw new IllegalArgumentException(
                    "This endposition is impossible. Try: startposition <= endposition < sequence.length");
        } else if (endpos - startpos + 1 != length) {
            throw new IllegalArgumentException("This sequence has not length " + length + ".");
        }
    }

    public double getFreq(int index, int observation) {
        double[] marginal = new double[4];
        for (int k = 0; k < condProb[index].length; k += 4) {
            for (int l = k; l < k + 4; l++) {
                marginal[l - k] += condProb[index][l];
            }
        }
        Normalisation.sumNormalisation(marginal);
        return marginal[observation];
    }

    /** @return true, if all positions in the given sequence from start to end (inclusive) are unobserved */
    public static boolean isCompletylUnObserved(Sequence<?> sequence, int start, int end) {
        for (int i = start; i <= end; i++) {
            if (sequence.discreteVal(i) != 4) {
                return false;
            }
        }
        return true;
    }

    /** Transforms this AlignmentBasedModel to a PhyloBayesModel */
    public PhyloBayesModel toPhyloBayesModel(String newickString) {
        PhyloBayesModel model = new PhyloBayesModel(this.length, this.order, newickString, this.alphabets);
        model.setCondProbs(this.getCondProbs());
        model.reinitParameters();
        return model;
    }
}
