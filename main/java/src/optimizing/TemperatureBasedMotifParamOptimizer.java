package optimizing;

import java.util.Arrays;

import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.LimitedMedianStartDistance;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.utils.SafeOutputStream;
import evolution.EvolModel;
import io.Alphabet;
import models.AbstractPhyloModel;
import util.MatrixLinearisation;

/**
 * This Optimizer is a prototype to learn <b>only</b> a PWM and if possible a filtering parameter from the data. It
 * assumes a simple star topology and the first species as reference species.
 * 
 * @author Martin Nettling
 * 
 */
public class TemperatureBasedMotifParamOptimizer extends AbstractFreeEnergyOptimizer {
    /** Mapping of parameters from position in motif to position in vectors */
    private int[] bufferPositions;
    /** Mapping of number of parameters from position in motif */
    private int[] bufferLength;

    private static double MIN_TEMPERATURE = 0.0;
    /**
     * Constructor.
     * 
     * @throws CloneNotSupportedException
     */
    public TemperatureBasedMotifParamOptimizer() throws CloneNotSupportedException {
        super(1e-4); // only possible with lambda-parametrisation
        tc = PreparedConditions.TC_MOTIF_ALL_POSITIONS.clone();
    }

    @Override
    public void setModel(AbstractPhyloModel model) {
        super.setModel(model);
        // check if a star toplogy was given
        if (bnh.getVirtualTree(0).numberOfNodes != bnh.getVirtualTree(0).numberOfLeafs + 1) {
            throw new IllegalArgumentException(
                    "A star topology is expected. In a star topology there exists only one non leaf node, the root");
        }
        
        this.dimension = 0; // one parameter is needed for the temperature
        this.bufferLength = new int[bnh.motifLength];
        this.bufferPositions = new int[bnh.motifLength];

        for (int i = 0; i < bnh.motifLength; i++) {
            bufferPositions[i] = this.dimension;
            bufferLength[i] = MatrixLinearisation.linearize(bnh.getVirtualTree(i).getEvolModel().getStatDistr()).length;
            dimension += bufferLength[i];
        }
        dimension += 1; // temperature
    }

    @Override
    public double evaluateFunction(double[] lambda) throws DimensionException, EvaluationException {
        initParameters (lambda);
        return getWeightedLogLikelihoodSum();
    }

    // Parameters needed in evaluateFunction

    
    private void initParameters(double[] lambda) {
        for (int i = 0; i < bnh.motifLength; i++) {
            EvolModel evol = bnh.getVirtualTree(i).getEvolModel();
            bnh.getVirtualTree(i).TEMPERATURE = Math.exp(lambda[lambda.length - 1]) + MIN_TEMPERATURE;
            
            double[] pi = new double[this.bufferLength[i]];
            // transform parameters to be usable in PhyloModels
            MatrixLinearisation.lambda2pi(
                    Arrays.copyOfRange(lambda, bufferPositions[i], bufferPositions[i] + bufferLength[i]), pi,
                    Alphabet.size);
            MatrixLinearisation.fillMatrix(pi, evol.getStatDistr()); // fills stationary distribution form vector pi
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
                    1e-3, new LimitedMedianStartDistance(10, .1), sus // output domain
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
        double[] lambda = new double[dimension];
        lambda[lambda.length - 1] = Math.log(bnh.getVirtualTree(0).TEMPERATURE - MIN_TEMPERATURE);
        int bufferPos = 0;

        for (int i = 0; i < bnh.motifLength; i++) {
            double[] tmp = MatrixLinearisation.pi2lambda(MatrixLinearisation.linearize(bnh.getVirtualTree(i)
                    .getEvolModel().getStatDistr()));
            for (int j = 0; j < tmp.length; j++) {
                lambda[bufferPos + j] = tmp[j];
            }
            bufferPos += tmp.length;
        }
        return lambda;
    }
}
