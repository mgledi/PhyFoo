package optimizing.branchlengths;

import java.util.ArrayList;

import algorithm.MeanFieldForBayesNet;
import bayesNet.BayesNetNode;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.LimitedMedianStartDistance;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.utils.SafeOutputStream;
import optimizing.AbstractFreeEnergyOptimizer;
import optimizing.PreparedConditions;

/**
 * This Optimizer is able to to optimize the branch lengths locally for a given set of {@link MeanFieldForBayesNet}s.
 * That means for each position in the underlying {@link MeanFieldForBayesNet}s the branches can be learned separately.
 * 
 * @author Chaos
 */
public class EdgeOptimizer extends AbstractFreeEnergyOptimizer {
    public static double MINIMUM_BRANCH_LENGTH = 0.0005;
    private boolean optimzeGlobally = false;
    private boolean optimizeBranchlengthsEqual = false;
    
    /** The position which should be trained */
    protected int pos;
    
    /**
     * Instantiates a new EdgeOptimizer. This Optimizer is able to train only one position of each element in the
     * underlying {@link MeanFieldForBayesNet}s at once.
     * 
     * @throws CloneNotSupportedException
     */
    public EdgeOptimizer(boolean optimzeGlobally, boolean equalBranchlengths) throws CloneNotSupportedException {
        super(1e-4); // only possible with lambda-parametrisation
        tc = PreparedConditions.SMALL_DIFFERENCE_OF_FUNCTIONS.clone();
        this.optimzeGlobally = optimzeGlobally;
        this.optimizeBranchlengthsEqual = equalBranchlengths;
    }

    @Override
    public double evaluateFunction(double[] edgesExp) throws DimensionException, EvaluationException {
        double[] edgesPi = new double[edgesExp.length];
        for (int i = 0; i < edgesPi.length; i++) {
            edgesPi[i] = (Math.exp(edgesExp[i]) + MINIMUM_BRANCH_LENGTH) / (1 + Math.exp(edgesExp[i]));
            if (Double.isNaN(edgesPi[i])) {
                edgesPi[i] = 1;
            }
        }
        // estimate positions to train on
        int start = 0, length = bnh.motifLength;
        if (!optimzeGlobally) {
            start = pos;
            length = pos + 1;
        }
        
        for (int i = start; i < length; i++) {
            int k = 0;
            for (int n = 0; n < bnh.getVirtualTree(pos).numberOfNodes; n++) {
                BayesNetNode bn = bnh.getVirtualTree(pos).getNode(n);
                if (!bn.props.isPhyloRoot()) {
                    bn.setDistanceToParent(edgesPi[k]);
                    // if branch length should be optimized independently, try to set next
                    if(!optimizeBranchlengthsEqual) {
                        k++;
                    }
                }
            }
            bnh.getVirtualTree(i).initParametersFromGF();
        }
        return this.getWeightedLogLikelihoodSum();
    }

    /** sets the position, which should be optimized */
    public void initPosition(int pos) {
        this.pos = pos;
    }

    /** startet die optimierung, ohne Ausgabe */
    public void startOptimizing() {
        this.startOptimizing(false);
    }

    public void startOptimizing(boolean output) {
        SafeOutputStream sus = getOutputStream(output);
        try {
            // independent of optimizeGlobally
            double avgPi = 0, count = 0;
            ArrayList<Double> edgesPi = new ArrayList<Double>();
            for (int i = 0; i < bnh.getVirtualTree(pos).numberOfNodes; i++) {
                BayesNetNode bn = bnh.getVirtualTree(pos).getNode(i);
                if (!bn.props.isPhyloRoot()) {
                    if(!optimizeBranchlengthsEqual) edgesPi.add(bn.getDistanceToParent());
                    avgPi += (bn.getDistanceToParent());
                    count ++;
                }
            }
            avgPi /= count;
            if(optimizeBranchlengthsEqual) {
                edgesPi.add(avgPi);
            }
            double[] edgesExp = new double[edgesPi.size()];
            dimension = edgesExp.length;
            for (int i = 0; i < edgesExp.length; i++) {
                double pi = edgesPi.get(i);
                if (pi == 1) {
                    pi -= 1e-6; // if edgelength is set to 1, I will get an numerical error
                }
                edgesExp[i] = Math.log(pi / (1 - pi));
            }
            sus.write("Train edges for position " + pos + "\n");
            Optimizer.optimize(Optimizer.QUASI_NEWTON_BFGS, new NegativeDifferentiableFunction(this), // initialize
                                                                                                      // function
                    edgesExp, // starting point of optimization
                    tc, // termination condition
                    1e-3, new LimitedMedianStartDistance(10, .1), sus // output domain
                    );

            sus.write("End Train edges for position " + pos + "\n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
