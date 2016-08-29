package optimizing;

import java.util.ArrayList;

import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.utils.Time;

public class IgnoringNegativeDifferenceOfFunctionCondition extends SmallDifferenceOfFunctionEvaluationsCondition {
    public static int ALLOWED_GLOBAL_NEGATIVE_DIFFERENCES = 5;
    public static int ALLOWED_CONSECUTIVE_NEGATIVE_DIFFERENCES = 10;
    private ArrayList<Double> differences;
   
    private double eps;
    
    public IgnoringNegativeDifferenceOfFunctionCondition(double epsilon) throws Exception {
        super(epsilon);
        differences = new ArrayList<Double>(250);
    }
    
    @Override
    public boolean doNextIteration( int iteration, double f_last, double f_current, double[] gradient, double[] direction, double alpha,
            Time t ) {
        differences.add(f_last - f_current);
        if(isGlobalConditionFulfilled() || isLocalConditionFulfilled()) {
            return false;
        } 
        return eps <= Math.abs(f_last-f_current);
    }
    
    private boolean isGlobalConditionFulfilled() {
        int count = 0;
        for(double d : differences) {
            if (d < 0) {
                count++;
            }
            if( count > ALLOWED_GLOBAL_NEGATIVE_DIFFERENCES) {
                return true;
            }
        }
        return false;
    }
    
    private boolean isLocalConditionFulfilled() {
        int count = 0;
        for(double d : differences) {
            if (d < 0) {
                count ++;
            } else {
                count = 0;
            }
            
            if( count > ALLOWED_CONSECUTIVE_NEGATIVE_DIFFERENCES) {
                return true;
            }
        }
        return false;
    }
    
    protected void set() {
        this.eps = (Double) parameter.getParameterAt(0).getValue();
    }
}
