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

public class MotifParamOptimizer extends AbstractFreeEnergyOptimizer {
    /** Mapping of parameters from position in motif to position in vectors */
    private int[] bufferPositions;
    /** Mapping of number of parameters from position in motif */
    private int[] bufferLength;

    /**
     * Constructor.
     * @throws CloneNotSupportedException 
     */
    public MotifParamOptimizer() throws CloneNotSupportedException {
        super(1e-4); // only possible with lambda-parametrisation
            tc = PreparedConditions.TC_MOTIF_ALL_POSITIONS.clone();
    }

    @Override
    public void setModel(AbstractPhyloModel model) {
        super.setModel(model);
        this.dimension = 0;
        this.bufferLength = new int[bnh.motifLength];
        this.bufferPositions = new int[bnh.motifLength];

        for (int i = 0; i < bnh.motifLength; i++) {
            bufferPositions[i] = this.dimension;
            bufferLength[i] = MatrixLinearisation.linearize(bnh.getVirtualTree(i).getEvolModel().getStatDistr()).length;
            dimension += bufferLength[i];
        }
    }

    @Override
    public double evaluateFunction(double[] lambda) throws DimensionException, EvaluationException {
        this.initParameters(lambda);
        return getWeightedLogLikelihoodSum();
    }

    private void initParameters(double[] lambda) {
        for (int i = 0; i < bnh.motifLength; i++) {
            EvolModel evol = bnh.getVirtualTree(i).getEvolModel();
            double[] tmp = new double[this.bufferLength[i]];
            MatrixLinearisation.lambda2pi(
                    Arrays.copyOfRange(lambda, bufferPositions[i], bufferPositions[i] + bufferLength[i]), tmp,
                    Alphabet.size);
            MatrixLinearisation.fillMatrix(tmp, evol.getStatDistr());
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
            Optimizer.optimize(
                    Optimizer.QUASI_NEWTON_BFGS,
                    new NegativeDifferentiableFunction(this), // initialize function
                    lambda, // starting point of optimization
                    tc, // termination condition
                    1e-3,
                    new LimitedMedianStartDistance(10, .1),
                    sus // output domain
                    );
            this.initParameters(lambda);
            sus.write("End Train \n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Returns the parameter-vector from the underlying BayesNet. The Vector
     * contains all Motif parameters
     */
    private double[] getParamLambdaVector() {
        double[] lambda = new double[dimension];
        int bufferPos = 0;

        for (int i = 0; i < bnh.motifLength; i++) {
            double[] tmp = MatrixLinearisation.pi2lambda(MatrixLinearisation.linearize(bnh.getVirtualTree(i).getEvolModel().getStatDistr()));
            for (int j = 0; j < tmp.length; j++) {
                lambda[bufferPos + j] = tmp[j];
            }
            bufferPos += tmp.length;
        }
        return lambda;
    }
}
