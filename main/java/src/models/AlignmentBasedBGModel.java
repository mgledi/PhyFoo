package models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

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

/**
 * This Model is the optimized version of {@link PhyloBackground} with branch lengths equal to 1.
 * 
 * @author Martin Nettling
 */
public class AlignmentBasedBGModel extends AbstractAlignmentBasedModel {

	private static Logger LOGGER = Logger.getLogger(PhyloBayesModel.class);

    public static boolean IGNORE_GAPS_IN_TRAINING = false;
    public static boolean IGNORE_GAPS_IN_LIKELIHOOD_CALC = false;
    public static double PSEUDO_COUNT = 1e-2;
    
    public AlignmentBasedBGModel(byte order, AlphabetContainer alphabet) {
        super(alphabet, order, 0);
        initPows();
        lnCondProb = new double[order + 1][];
        condProb = new double[order + 1][];

        for (int m = 0; m <= order; m++) {
            lnCondProb[m] = new double[pows[m] * 4];
            condProb[m] = new double[pows[m] * 4];
        }
        fillLnCondProbRandomly();
    }

    public AlignmentBasedBGModel(StringBuffer xml) throws NonParsableException {
        super(xml);
        initPows();
    }

    @Override
    public void train(Sample data, double[] weights) throws Exception {
        double[/* motif-position */][/* base pair-hash */] counts = new double[order + 1][];
        for (int i = 0; i <= order; i++) {
            counts[i] = new double[lnCondProb[i].length];
            Arrays.fill(counts[i], PSEUDO_COUNT);
        }

        Sequence<?>[] sequences = data.getAllElements();
        for (int w = 0; w < sequences.length; w++) {
            MultiDimensionalDiscreteSequence parent;
            int start = 0, end = 0;

            if (sequences[w] instanceof MultiDimensionalDiscreteSequence) {
                parent = (MultiDimensionalDiscreteSequence) sequences[w];
                end = parent.getLength();
            } else {
                RecursiveSequence<?> subseq = (RecursiveSequence<?>) sequences[w];
                parent = (MultiDimensionalDiscreteSequence) subseq.getParent();
                start = subseq.getIndex(0);
                subseq.getIndex(subseq.getLength());
            }
            // run over all positions in the given sequence (length)
            ArrayList<int[]> observations = new ArrayList<int[]>();
            for (int o = 0; o < parent.getNumberOfSequences(); o++) {
                for (int l = start; l < end; l++) {
                    // k,o is the position in the sequence -> determine count hash
                    observations.clear();
                    int preSize = 0;
                    for (int k = 0; k <= order && k <= (l - start); k++) {
                        if (parent.containerForPhyloBayes[l - k][o] == 4 && IGNORE_GAPS_IN_TRAINING) {
                            observations.clear();
                            break;
                        } else if (parent.containerForPhyloBayes[l - k][o] == 4) {
                            int oldPreSize = preSize;
                            preSize = observations.size();
                            if (observations.size() == 0) {
                                for (int a = 0; a < 4; a++) {
                                    int[] tmp = new int[order + 1];
                                    tmp[k] = a;
                                    observations.add(tmp);
                                }
                            } else {
                                for (int i = oldPreSize; i < preSize; i++) {
                                    for (int a = 0; a < 4; a++) {
                                        int[] tmp = Arrays.copyOf(observations.get(i), order + 1);
                                        tmp[k] = a;
                                        observations.add(tmp);
                                    }
                                }
                            }
                        } else {
                            if (observations.size() == 0) {
                                int[] tmp = new int[order + 1];
                                observations.add(tmp);
                            }
                            for (int[] obs : observations) {
                                obs[k] = parent.containerForPhyloBayes[l - k][o];
                            }
                        }
                    }
                    // increment count-hash by the weight[w]
                    int size = observations.size() - preSize;
                    for (; preSize < observations.size(); preSize ++) {
                        int[] obs = observations.get(preSize);
                        int hash = 0;
                        for (int k = 0; k <= order && k <= (l - start); k++) {
                            hash += pows[k] * obs[k];
                        }
                        // if a gap was observed the weight must been distributed between all observations
                        counts[Math.min(l, order)][hash] += (weights == null ? 1 : weights[w]) / size;
                    }
                }
            }
        }
        // Normalisation
        for (int i = 0; i <= order; i++) {
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
     * Sets the conditional probability for position i
     */
    @Override
    public void setLnCondProb(double[] lambda, int position) {
        if (lambda.length != lnCondProb[position].length) {
            throw new IllegalArgumentException("The given lambda array must be of length "
                    + lnCondProb[position].length + ".");
        } else if (position > length) {
            throw new IllegalArgumentException("The given position must be smaller than the length (" + length
                    + "), but is " + position + ".");
        } else {
            lnCondProb[position] = Arrays.copyOf(lambda, lambda.length);
        }
    }

    /** @return a pointer to the underlying conditional probability array at position i */
    @Override
    public double[] getLnCondProb(int i) {
        return lnCondProb[i];
    }

    public double[] getCondProb(int i) {
        return condProb[i];
    }

    /**
     * Estimates the conditional probability for nucleotide at position posSeq in the given sequence. The
     * log-probability.
     * 
     * @return
     */
    public double getLnCondProb(Sequence<?> sequence, final int posSeq) {
        if (AlignmentBasedModel.isCompletylUnObserved(sequence, Math.max(0, posSeq - order), posSeq)) {
            return Math.log(0.25);
        }
        // prepare all int[] - arrays, which are represented by the subsequence at posSeq
        ArrayList<int[]> observations = new ArrayList<int[]>();
        int preSize = 0;
        for (int k = 0; k <= order && k <= posSeq; k++) {
            if (sequence.discreteVal(posSeq - k) == 4 && IGNORE_GAPS_IN_LIKELIHOOD_CALC) {
                return 0; // HACKY
            } else if (sequence.discreteVal(posSeq - k) == 4) { // if gap is observed, create all possible observations 
                preSize = observations.size();
                if (observations.size() ==  0) { // if no observation was added yet
                    for (int a = 0; a < 4; a++) {
                        int[] tmp = new int[order + 1];
                        tmp[k] = a;
                        observations.add(tmp);
                    }
                } else { // if more than one observation is available
                    // copy all (not yet complete) observations and add any base
                    for (int i = 0; i < preSize; i++) {
                        for (int a = 0; a < 4; a++) {
                            int[] tmp = Arrays.copyOf(observations.get(i), order + 1);
                            tmp[k] = a;
                            observations.add(tmp);
                        }
                    }
                }
            } else { // if concrete symbol is observed, add symbol to all observations
                if (observations.size() == 0) {
                    int[] tmp = new int[order + 1];
                    observations.add(tmp);
                }

                for (int[] obs : observations) {
                    obs[k] = sequence.discreteVal(posSeq - k);
                }
            }
        }
        // run over all possible observations and add up probabilities
        double sumProb = 0; 
        int size = observations.size() - preSize;
        for (; preSize < observations.size(); preSize ++) {
            int[] obs = observations.get(preSize);
            int access = 0;
            for (int k = 0; k <= order && k <= posSeq; k++) {
                access += pows[k] * obs[k];
            }
            sumProb += condProb[Math.min(posSeq, order)][access] / size;
        }
        return Math.log(sumProb);
    }

    /* (non-Javadoc)
     * 
     * @see de.jstacs.models.AbstractModel#getLogProbFor(de.jstacs.data.Sequence, int, int) */
    @Override
    public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws NotTrainedException, Exception {
        if (endpos < startpos) {
            return 0;
        }
        check(sequence, startpos, endpos);

        MultiDimensionalDiscreteSequence alignment = (MultiDimensionalDiscreteSequence) sequence;
        double logSum = 0;
        for (int l = startpos; l <= endpos; l++) {
            for (int o = 0; o < alignment.getNumberOfSequences(); o++) {
                logSum += this.getLnCondProb(alignment.getSequence(o), l);
            }
        }

        // sum up all nucleotides and probabilities
        return logSum;
    }

    @Override
    public AlignmentBasedBGModel clone() {
        StringBuffer xml = this.toXML();
        try {
            return new AlignmentBasedBGModel(xml);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return this;
    }

    @Override
    public String getInstanceName() {
        return AlignmentBasedBGModel.class.getSimpleName();
    }

    @Override
    public StringBuffer toXML() {
        StringBuffer xml = new StringBuffer();
        XMLParser.appendObjectWithTags(xml, order, "order");
        XMLParser.appendObjectWithTags(xml, alphabets, "alphabets");
        XMLParser.appendObjectWithTags(xml, lnCondProb, "lnCondProb");
        XMLParser.appendObjectWithTags(xml, condProb, "condProb");
        return xml;
    }

    @Override
    protected void fromXML(StringBuffer xml) throws NonParsableException {
        length = 0;
        order = (Byte) XMLParser.extractObjectForTags(xml, "order");
        alphabets = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "alphabets");
        lnCondProb = (double[][]) XMLParser.extractObjectForTags(xml, "lnCondProb");
        condProb = (double[][]) XMLParser.extractObjectForTags(xml, "condProb");
    }

    public byte getMaximalMarkovOrder() {
        return (byte) (order);
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
        }
    }

    // this caching speedups the Calculation of getLogProbFor by factor of five
    private HashMap<Sequence, HashMap<Integer, HashMap<Integer, Double>>> _cachedScores = new HashMap<Sequence, HashMap<Integer, HashMap<Integer, Double>>>();

    /** cleans the cache */
    public void cleanCache() {
        _cachedScores = new HashMap<Sequence, HashMap<Integer, HashMap<Integer, Double>>>();
    }

    /** true: if a score was already computed */
    public boolean isScoreCached(Sequence sequence, int startpos, int endpos) {
        if (_cachedScores.containsKey(sequence) && _cachedScores.get(sequence).containsKey(startpos)
                && _cachedScores.get(sequence).get(startpos).containsKey(endpos)) {
            return true;
        } else {
            return false;
        }
    }

    /** inits needed HashMaps for storing a score */
    public void prepareCache(Sequence sequence, int startpos, int endpos, double score) {
        if (!_cachedScores.containsKey(sequence)) {
            _cachedScores.put(sequence, new HashMap<Integer, HashMap<Integer, Double>>());
        }

        if (!_cachedScores.get(sequence).containsKey(startpos)) {
            _cachedScores.get(sequence).put(startpos, new HashMap<Integer, Double>());
        }

        _cachedScores.get(sequence).get(startpos).put(endpos, score);
    }

    /** Transforms this AlignmentBasedBGModel to a PhyloBackground 
     * @throws Exception */
    public PhyloBackground toPhyloBackground(String newickString) throws Exception {
        PhyloBackground model = new PhyloBackground(this.order, newickString, this.alphabets);
        model.setCondProbs(this.getCondProbs());
        model.reinitParameters();
        return model;
    }
}
