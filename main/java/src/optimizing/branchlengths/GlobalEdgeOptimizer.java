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
 * This Optimizer is able to to optimize the branch lengths globally for a given set of {@link MeanFieldForBayesNet}s.
 * 
 * @author Chaos
 */
public class GlobalEdgeOptimizer extends AbstractFreeEnergyOptimizer {
    
    /**
     * Constructor. Instantiates a new {@link GlobalEdgeOptimizer}.
     * @throws CloneNotSupportedException 
     */
    public GlobalEdgeOptimizer() throws CloneNotSupportedException  {
        super(1e-4); // only possible with lambda-parameterization
        tc = PreparedConditions.SMALL_DIFFERENCE_OF_FUNCTIONS.clone();
    }

    @Override
    public double evaluateFunction(double[] edgesExp) throws DimensionException, EvaluationException {
        double[] edgesPi = new double[edgesExp.length];
        for (int i = 0; i < edgesPi.length; i++) {
            edgesPi[i] = (Math.exp(edgesExp[i]) + EdgeOptimizer.MINIMUM_BRANCH_LENGTH) / (1 + Math.exp(edgesExp[i]));
            if (Double.isNaN(edgesPi[i])) {
                edgesPi[i] = 1;
            }
        }

        for (int i = 0; i < bnh.motifLength; i++) {
            int k = 0;
            for (int n = 0; n < bnh.getVirtualTree(i).numberOfNodes; n++) {
                BayesNetNode bn = bnh.getVirtualTree(i).getNode(n);
                if (!bn.props.isPhyloRoot()) {
                    bn.setDistanceToParent(edgesPi[k++]);
                }
            }
            bnh.getVirtualTree(i).initParametersFromGF();
        }
        double val = this.getWeightedLogLikelihoodSum();
//        double val = this.getSumKLDivergences();
//        double val = this.getNegativeSumDifferenceOfIC();
//        double val = this.getSumEuclideanDistance();
        return val;
    }

    /** Starts optimizing without output */
    public void startOptimizing() {
        this.startOptimizing(false);
    }

    /** Starts optimizing */
    public void startOptimizing(boolean output) {
        try {
            this.initSpeciesSamples(); //TODO: only init if needed
        } catch (Exception e) {
            throw new RuntimeException("Could not init single species samples",e);
        }
        
        SafeOutputStream sus = getOutputStream(output);
        try {
            // load all edges into array
            ArrayList<Double> edgesPi = new ArrayList<Double>();
            for (int i = 0; i < bnh.getVirtualTree(0).numberOfNodes; i++) {
                BayesNetNode bn = bnh.getVirtualTree(0).getNode(i);
                if (!bn.props.isPhyloRoot()) {
                    edgesPi.add(bn.getDistanceToParent());
                }
            }
            double[] edgesExp = new double[edgesPi.size()];
            dimension = edgesExp.length;
            for (int i = 0; i < edgesExp.length; i++) {
                double pi = edgesPi.get(i) == 1 ? edgesPi.get(i) - 1e-6 : edgesPi.get(i);
                edgesExp[i] = Math.log(pi / (1 - pi));
            }
            sus.write("Train edges globally \n");
            Optimizer.optimize(
                    Optimizer.QUASI_NEWTON_BFGS,
                    new NegativeDifferentiableFunction(this), // initialize function
                    edgesExp, // starting point of optimization
                    tc, // termination condition
                    1e-3,
                    new LimitedMedianStartDistance(10, .1),
                    sus // output domain
                    );

            sus.write("End Train edges globally \n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}