package optimizing;

import algorithm.MeanFieldForBayesNet;
import de.jstacs.algorithms.optimization.DimensionException;
import de.jstacs.algorithms.optimization.EvaluationException;
import de.jstacs.algorithms.optimization.LimitedMedianStartDistance;
import de.jstacs.algorithms.optimization.NegativeDifferentiableFunction;
import de.jstacs.algorithms.optimization.Optimizer;
import de.jstacs.utils.SafeOutputStream;
import evolution.EvolModel;
import io.Alphabet;
import util.MatrixLinearisation;

/**
 * An optimizer for free Energy of a BayesNet, concerning the equilibrium frequencies of the evolutionary models in the
 * underlying {@link MeanFieldForBayesNet}s. This optimizer handles each position separately.
 */
public class FE_byParameter_ForBayesNodeJstacs extends AbstractFreeEnergyOptimizer {
    /** The position which should be trained */
    private int pos;

    /**
     * Instantiates a new Optimizer.
     * @throws CloneNotSupportedException 
     */
    public FE_byParameter_ForBayesNodeJstacs() throws CloneNotSupportedException {
        super(1e-4); // only possible with lambda-parametrisation
        tc = PreparedConditions.TC_MOTIF_SEPARATED_POSITIONS.clone();
    }


    /** sets the position, which should be optimized */
    public void initPosition(int pos) {
        this.pos = pos;
    }

    private double[] tmpPi;

    @Override
    /**
     * berechnet einen Funktionswert zu dem �bergebenen Vector
     */
    public double evaluateFunction(double[] lambda) throws DimensionException, EvaluationException {
        initParameters(lambda);
        return getWeightedLogLikelihoodSum();
    }

    private void initParameters(double[] lambda) {
        EvolModel evol = bnh.getVirtualTree(pos).getEvolModel(); // get a pointer to the FelsteinModel
        MatrixLinearisation.lambda2pi(lambda, tmpPi, Alphabet.size);
        MatrixLinearisation.fillMatrix(tmpPi, evol.getStatDistr());
        bnh.getVirtualTree(pos).initParametersFromGF();
    }
    
    /** startet die optimierung, ohne Ausgabe */
    public void startOptimizing() {
        this.startOptimizing(false);
    }
    
    /** startet die optimierung */
    public void startOptimizing(boolean output) {
        SafeOutputStream sus = getOutputStream(output);
        try {
            double[] tmp = MatrixLinearisation.linearize(bnh.getVirtualTree(pos).getEvolModel().getStatDistr());
            this.dimension = tmp.length;
            double[] lambda = MatrixLinearisation.pi2lambda(tmp);
            tmpPi = new double[lambda.length];
            sus.write("Train position " + pos + "\n");
            Optimizer.optimize(
                    Optimizer.QUASI_NEWTON_BFGS,
                    new NegativeDifferentiableFunction(this), // initialize function
                    lambda, // starting point of optimization
                    tc, // termination condition
                    1e-3,
                    new LimitedMedianStartDistance(10, .1),
                    sus // output domain
                    );
            initParameters(lambda);
            sus.write("End Train position " + pos + "\n");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
