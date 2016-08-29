package models;

import java.util.Arrays;
import java.util.HashMap;

import algorithm.MeanFieldForBayesNet;
import bayesNet.BayesNet;
import bayesNet.BayesNetHandler;
import bayesNet.CPF;
import bayesNet.VirtualTree;
import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import evolution.EvolModel;
import io.Alphabet;
import util.MatrixLinearisation;
import util.Util;

public abstract class AbstractPhyloModel extends PhyloPreparedAbstractModel {
    /**
     * false = Each branch is learned independently
     * true = all branch lengths will have the same length
     */
    public static boolean EDGE_LEARNING_ONE_LENGTH = false;

    /**
     * indicates to cache the calculated MeanFieldForBayesNet-Objects. This brings performance much performance but it
     * needs also much memory. So enable this for "small" Training-sets
     */
    public static boolean ENABLE_CACHING = true;

    /** a handler for the underlying bayesnet */
    protected BayesNetHandler bnh;

    transient private HashMap<Sequence<?>, MeanFieldForBayesNet> mfCache = new HashMap<Sequence<?>, MeanFieldForBayesNet>();

    /**
     * @param alphabets
     * @param order
     * @param length
     */
    public AbstractPhyloModel(AlphabetContainer alphabets, byte order, int length) {
        super(alphabets, order, length);
    }

    /**
     * @param xml
     * @throws NonParsableException
     */
    public AbstractPhyloModel(StringBuffer xml) throws NonParsableException {
        super(xml);
    }

    /*
     * (non-Javadoc)
     * 
     * @see de.jstacs.models.Model#train(de.jstacs.data.Sample)
     */
    public void train(Sample data) throws Exception {
        double[] weights = new double[data.getAllElements().length];
        Arrays.fill(weights, 1);
        train(data, weights);
    }

    /** @return Pointer to the used {@link BayesNet} */
    public BayesNetHandler getBNH() {
        return bnh;
    }

    /** reinitialise all parameters. Normally a call of this method is not needed */
    public void reinitParameters() {
        for (int i = 0; i < bnh.motifLength; i++) {
            bnh.getVirtualTree(i).initParametersFromGF();
        }
    }

    /** Reinitializes for each position in the branch lengths of the underlying phylogenetic trees. */
    public void reinitEdgelengths(String... trees) {
        if (trees.length != bnh.getVirtualTrees().length && trees.length != 1) {
            throw new IllegalArgumentException("The number of given trees must be equal to length of this model");
        }
        if (trees.length == 1) {
            for (int i = 0; i <= (this.getLength() == 0 ? order : this.getLength() - 1); i++) {
                bnh.getVirtualTree(i).reinitEdglengths(trees[0]);
            }
        } else {
            for (int i = 0; i < trees.length; i++) {
                bnh.getVirtualTree(i).reinitEdglengths(trees[i]);
            }
        }
    }

    /**
     * handles the cache for all Meanfieldcalculations
     * 
     * @throws Exception
     */
    @SuppressWarnings({ "unchecked", "rawtypes" })
    protected MeanFieldForBayesNet getMFfromCache(Sequence sequence) {
        MeanFieldForBayesNet tmpPointer;
        // if caching is not enabled, always instantiate new MeanField
        if (!ENABLE_CACHING) {
            tmpPointer = new MeanFieldForBayesNet(bnh, sequence);
            tmpPointer.initObservation();
            tmpPointer.init();
            return tmpPointer;
        } else {
            if (!mfCache.containsKey(sequence.hashCode())) {
                tmpPointer = new MeanFieldForBayesNet(bnh, sequence);
                tmpPointer.initObservation();
                tmpPointer.init();
                mfCache.put(sequence, tmpPointer);
            }
            tmpPointer = mfCache.get(sequence);
            tmpPointer.initObservation();
            return tmpPointer;
        }
    }

    @Override
    public void setLnCondProb(double[] lambda, int position) {
        double[] condProb = new double[lambda.length];
        for (int a = 0; a < condProb.length; a++) {
            condProb[a] = Math.exp(lambda[a]);
        }
        this.setCondProb(condProb, position);
    }

    @Override
    public void setCondProb(double[] condProb, int position) {
        MatrixLinearisation.fillMatrix(condProb, getBNH().getVirtualTree(position).getEvolModel().getStatDistr());
        getBNH().getVirtualTree(position).initParametersFromGF();
    }

    /** @return a copy of the underlying logarithmic conditional probability array at position i */
    @Override
    public double[] getLnCondProb(int i) {
        double[] lnCondProb = getCondProb(i);
        for (int a = 0; a < lnCondProb.length; a++) {
            lnCondProb[a] = Math.log(lnCondProb[a]);
        }
        return lnCondProb;
    }

    /** @return a copy of the underlying logarithmic conditional probabilities */
    @Override
    public double[][] getLnCondProbs() {
        double[][] lnCondProbs = getCondProbs();
        for (int i = 0; i < lnCondProbs.length; i++) {
            for (int a = 0; a < lnCondProbs[i].length; a++) {
                lnCondProbs[i][a] = Math.log(lnCondProbs[i][a]);
            }
        }
        return lnCondProbs;
    }

    /** @return a copy of the underlying conditional probability array at position i */
    @Override
    public double[] getCondProb(int i) {
        return MatrixLinearisation.linearize(bnh.getVirtualTree(i).getEvolModel().getStatDistr());

    }

    /** @return a copy of the underlying conditional probability array for species s */
    @Override
    public double[][] getCondProbsPerSpecies(int s) {
        if (s < 0 && s >= bnh.getVirtualTree(0).numberOfLeafs) {
            throw new IllegalArgumentException("The given index " + s
                    + " does not correspond to a species in the underlying tree.");
        }
        if(order > 0) {
            throw new UnsupportedOperationException("This method is not implemented for order > 0");
        }
        if(bnh.getVirtualTree(0).numberOfLeafs != bnh.getVirtualTree(0).numberOfNodes - 1) {
            throw new UnsupportedOperationException("This method is not implemented for general topologies.");
        }
        double[][] m = new double[this.length == 0 ? order + 1 : this.length][];
        boolean rememberMe = CPF.ALLOW_RETURN_INDEPENDENT_TRANSITION;
        CPF.ALLOW_RETURN_INDEPENDENT_TRANSITION = false;
        
        this.reinitParameters();
        for (int i = 0; i < m.length; i++) {
            double[] root_pi = bnh.getVirtualTree(i).getNode(0).CPF.getCondProb()[0];
            double[][] FS81 = bnh.getVirtualTree(i).getLeaf(s).CPF.getCondProb();
//            Util.visualizeMatrix(FS81);
            m[i] = VirtualTree.times(root_pi, FS81);
        }
        
        CPF.ALLOW_RETURN_INDEPENDENT_TRANSITION = rememberMe;
        return m;
    }

    /** @return a copy of the underlying conditional probability array */
    @Override
    public double[][] getCondProbs() {
        double[][] m = new double[this.length == 0 ? order + 1 : this.length][];
        for (int i = 0; i < m.length; i++) {
            m[i] = getCondProb(i);
        }
        return m;
    }

    @Override
    public void fillLnCondProbRandomly() {
        EvolModel gfs = EvolModel.getInstance(1);
        for (int i = 0; i < bnh.getVirtualTrees().length; i++) {
            gfs = EvolModel.getInstance((int) Math.pow(Alphabet.size, bnh.getVirtualTree(i).getNode(0).numberOfParents));
            bnh.getVirtualTree(i).setEvolModel(gfs);

            bnh.getVirtualTree(i).getEvolModel()
                    .setStatDistr(Util.getRandomStochMatrix(gfs.getDimension(), Alphabet.size));
            bnh.getVirtualTree(i).initParametersFromGF();
        }
    }
}
