package optimizing;

import models.AbstractPhyloModel;
import bayesNet.BayesNetHandler;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.LimitedMedianStartDistance;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.utils.SafeOutputStream;

/**
 * This Optimizer is a prototype to learn <b>only</b> a PWM and if possible a filtering parameter from the data. It
 * assumes a simple star topology and the first species as reference species.
 * 
 * @author Martin Nettling
 * 
 */
public class TemperatureOptimizer extends AbstractFreeEnergyOptimizer {
    private static double MIN_TEMPERATURE = 0.0;
    private boolean learnLocalTemperatures = false;    
    
    /**
     * Constructor.
     * @param localTemperaturs 
     * 
     * @throws CloneNotSupportedException
     */
    public TemperatureOptimizer(boolean learnLocalTemperatures) throws CloneNotSupportedException {
        super(1e-4); // only possible with lambda-parametrisation
        tc = PreparedConditions.TC_MOTIF_ALL_POSITIONS.clone();
        this.learnLocalTemperatures = learnLocalTemperatures;
    }

    @Override
    public void setModel(AbstractPhyloModel model) {
        super.setModel(model);
        // check if a star toplogy was given
        if (bnh.getVirtualTree(0).numberOfNodes != bnh.getVirtualTree(0).numberOfLeafs + 1) {
            throw new IllegalArgumentException(
                    "A star topology is expected. In a star topology there exists only one non leaf node, the root");
        }
        if(learnLocalTemperatures) {
            this.dimension = bnh.motifLength; // one parameter is needed for the temperature
        } else {
            this.dimension = 1; // one parameter is needed for the temperature
        }
    }

    @Override
    public double evaluateFunction(double[] lambda) throws DimensionException, EvaluationException {
        initParameters (lambda);
        return getWeightedLogLikelihoodSum();
    }

    // Parameters needed in evaluateFunction

    
    private void initParameters(double[] lambda) {
        for (int i = 0; i < bnh.motifLength; i++) {
            if(learnLocalTemperatures) {
                bnh.getVirtualTree(i).TEMPERATURE = Math.exp(lambda[i]) + MIN_TEMPERATURE;
            } else {
                bnh.getVirtualTree(i).TEMPERATURE = Math.exp(lambda[0]) + MIN_TEMPERATURE;
            }
            bnh.getVirtualTree(i).initParametersFromGF();
        }
    }
    
    public void startOptimizing() {
        this.startOptimizing(false);
    }

    public void startOptimizing(boolean output) {
        SafeOutputStream sus = getOutputStream(output);
        try {
            double[] lambda = this.getParamLambdaVector();

            sus.write("Train all Positions \n");
            Optimizer.optimize(Optimizer.QUASI_NEWTON_BFGS, new NegativeDifferentiableFunction(this), // initialize
                                                                                                      // function
                    lambda, // starting point of optimization
                    tc, // termination condition
                    1e-3, new LimitedMedianStartDistance(10, .01), sus // output domain
                    );
            initParameters(lambda);
            sus.write("End Train \n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Returns the parameter-vector from the underlying BayesNet. The Vector contains all Motif parameters
     */
    private double[] getParamLambdaVector() {
        if(learnLocalTemperatures) {
            double[] lambda = new double[bnh.motifLength];
            for (int i = 0; i < lambda.length; i++) {
                lambda[i] = Math.log(bnh.getVirtualTree(0).TEMPERATURE - MIN_TEMPERATURE);
            }
            return lambda;
        } else {
            return new double[]{Math.log(bnh.getVirtualTree(0).TEMPERATURE - MIN_TEMPERATURE)};
        }
    }
}

