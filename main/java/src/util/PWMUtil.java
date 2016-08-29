package util;

import java.util.Random;

import bayesNet.BayesNetHandler;
import de.jstacs.models.AbstractModel;
import de.jstacs.models.Model;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.discrete.inhomogeneous.InhCondProb;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.utils.Normalisation;
import models.AlignmentBasedModel;
import models.PhyloPreparedAbstractModel;

public class PWMUtil {
    /** calculates the reverse pwm from a given pwm */
    public static double[][] getReversePWM(double[][] pwm) {
        double[][] reversePWM = new double[pwm.length][4];
        int length = pwm.length - 1;
        int alph = 3;
        for (int i = 0; i < pwm.length; i++) {
            for (int a = 0; a < 4; a++) {
                reversePWM[length - i][alph - a] = pwm[i][a];
            }
        }
        return reversePWM;
    }

    /**
     * calculates a pwm from marginal distributions in the root notes of a
     * bayesnet
     */
    public static double[][] getPWMFromBayesNet(BayesNetHandler bnh) {
        double[][] pwm = new double[bnh.motifLength][4];
        bnh.calculateCombinedProbs();
        for (int i = 0; i < bnh.motifLength; i++) {
            pwm[i] = bnh.getVirtualTree(i).getNode(0).CPF.getFullMarginalizedDistr();
        }
        return pwm;
    }

    /**
     * calculates a pwm from marginal distributions in the root notes of a
     * bayesnet
     */
    public static double[][] getPWMFrom(BayesianNetworkModel bnh) {
        InhCondProb[] constraints = bnh.constraints;
        double[][] pwm = new double[bnh.getLength()][4];
        for (int i = 0; i < constraints.length; i++) {
            for (int a = 0; a < 4; a++) {
                pwm[i][a] = constraints[i].getFreq(a);
            }
        }
        return pwm;
    }

    public static double[][] getPWMFrom(AlignmentBasedModel fg) {
        double[][] pwm = new double[fg.getLength()][4];
        for (int i = 0; i < fg.getLength(); i++) {
            for (int a = 0; a < 4; a++) {
                pwm[i][a] = fg.getFreq(i, a);
            }
        }
        return pwm;
    }

    public static double getMeanPolarity(double[][] PWM) {
        double avgPol = 0;
        for (int i = 0; i < PWM.length; i++) {
            double max = 0;
            for (int a = 0; a < PWM[i].length; a++) {
                if (PWM[i][a] > max) {
                    max = PWM[i][a];
                }
            }
            avgPol += max;
        }
        return avgPol / PWM.length;
    }

    
    /**
     * Parses a PWM from the given comma separated string
     * 
     * @param sPWM
     * @return the parsed PWM
     */
    public static double[][] parsePWMMirrored(String sPWM) {
        String[] explode = sPWM.split(",");
        double[][] pwm = new double[explode.length / 4][4];
        int l = 0;
        for (int a = 0; a < 4; a++) {
            for (int i = 0; i < pwm.length; i++) {
                    pwm[i][a] = Double.valueOf(explode[l]);
                    l++;
            }
        }
        for (int i = 0; i < pwm.length; i++) {
            Normalisation.sumNormalisation(pwm[i]);
        }
        return pwm;
    }
    
    /**
     * Parses a PWM from the given comma separated string
     * 
     * @param sPWM
     * @return the parsed PWM
     */
    public static double[][] parsePWM(String sPWM) {
        String[] explode = sPWM.split(",");
        double[][] pwm = new double[explode.length / 4][4];
        int l = 0;
        for (int i = 0; i < explode.length;) {
            for (int a = 0; a < 4; a++) {
                pwm[l][a] = Double.valueOf(explode[i]);
                i++;
            }
            Normalisation.sumNormalisation(pwm[l]);
            l++;
        }
        return pwm;
    }

    /**
     * This method returns an array of int that reflects a random order of the given length.
     * 
     * @param length
     * @return
     */
    public static int[] getRandomOrder(int length) {
        int[] order = new int[length];
        for (int i = 0; i < order.length; i++) {
            order[i] = i;
        }
        Random r = new Random();
        for (int i = 0; i < order.length / 2; i++) {
            int swap = r.nextInt(order.length - i) + i;
            int tmp = order[i];
            order[i] = order[swap];
            order[swap] = tmp;
        }
        return order;
    }
    
    public static void setPWM(BayesianNetworkModel model, double[][] pwm) {
        if(model.getLength() != pwm.length) {
            throw new IllegalArgumentException("Model and PWM must have the same length");
        }
        for (int i = 0; i < pwm.length; i++) {
            if (((BayesianNetworkModel) model).constraints[i].freq.length != 4) {
                throw new IllegalArgumentException("Can return PWM only for models with order 0");
            }
            for(int a=0; a < pwm[i].length; a++) {
                ((BayesianNetworkModel) model).constraints[i].freq[a] = pwm[i][a];
                ((BayesianNetworkModel) model).constraints[i].lnFreq[a] = Math.log(pwm[i][a]);
            }
        }
    }

	public static double[][] getPWM(Model model) {
		return getCondProbs(model, true);
	}
    /**
     * Returns the PWM from an {@link StrandModel}, {@link AbstractModel} or {@link BayesianNetworkModel}
     * 
     * @param model
     * @return
     */
	public static double[][] getCondProbs(Model model, boolean pwmOnly) {
        if (model instanceof StrandModel) {
            model = ((StrandModel) model).model[0];
        }

        if (model.getLength() == 0) {
			throw new IllegalArgumentException("Can return CondProbs only for models with length > 0");
        }

        if (model instanceof PhyloPreparedAbstractModel) {
			if (((PhyloPreparedAbstractModel) model).getOrder() > 0 && pwmOnly) {
                throw new IllegalArgumentException("Can return PWM only for models with order 0");
            }
            return Util.arraycopy(((PhyloPreparedAbstractModel) model).getCondProbs());
        } else if (model instanceof BayesianNetworkModel) {
            double[][] pwm = new double[model.getLength()][];
            for (int i = 0; i < pwm.length; i++) {
				if (((BayesianNetworkModel) model).constraints[i].freq.length != 4 && pwmOnly) {
                    throw new IllegalArgumentException("Can return PWM only for models with order 0");
                }
                pwm[i] = ((BayesianNetworkModel) model).constraints[i].freq;
            }
            return pwm;
        }
        throw new IllegalArgumentException("Unknown model " + model.getClass());
    }
}
