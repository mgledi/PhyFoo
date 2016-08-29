package optimizing;

import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.SmallGradientConditon;
import de.jstacs.algorithms.optimization.termination.SmallStepCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;

/**
 * This class contains several {@link TerminationCondition}s. They will be initialized in a static manner. These
 * conditions are used during optimizing by cloning them. This is much faster than reinstantiating.
 * 
 * @author Martin Nettling
 */
public class PreparedConditions {
    public static SmallDifferenceOfFunctionEvaluationsCondition SMALL_DIFFERENCE_OF_FUNCTIONS;
    public static CombinedCondition TC_MOTIF_SEPARATED_POSITIONS;
    public static CombinedCondition TC_MOTIF_ALL_POSITIONS;
    static {
        try {
            TC_MOTIF_SEPARATED_POSITIONS = new CombinedCondition(3,
                    new SmallDifferenceOfFunctionEvaluationsCondition(1e-4),
                    new IterationCondition(4), 
                    new SmallGradientConditon(1e-4), 
                    new SmallStepCondition(1e-4));
            TC_MOTIF_ALL_POSITIONS = new CombinedCondition(2, 
                    new SmallDifferenceOfFunctionEvaluationsCondition(1e-5),
                    new SmallGradientConditon(1e-4), 
                    new SmallStepCondition(1e-4));
            SMALL_DIFFERENCE_OF_FUNCTIONS = new SmallDifferenceOfFunctionEvaluationsCondition(1e-3);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
