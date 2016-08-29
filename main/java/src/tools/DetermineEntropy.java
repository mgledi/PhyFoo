package tools;

import java.io.File;
import java.util.Properties;

import de.jstacs.io.FileManager;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import models.AbstractPhyloModel;
import models.PhyloPreparedAbstractModel;
import util.Config;
import util.EvaluationUtil;

/**
 * Loads a model and extracts the underlying PWM. Then calculates the information content of this PWM and prints it on
 * stdout.
 * 
 * @author Martin Nettling
 * 
 */
public class DetermineEntropy {
    private static Properties PROPS = new Properties();

    public static void main(String[] args) throws Exception {
        Config.parseProperties(PROPS, args, false);
        String modelFile = Config.BASE_DIR + "/" + args[0];
        int species = Integer.valueOf(PROPS.getProperty("SPECIES", "-1"));
        StringBuffer xml = FileManager.readFile(new File(modelFile));
        PhyloPreparedAbstractModel model = ((PhyloPreparedAbstractModel) (new SingleHiddenMotifMixture(xml)).model[0]);
        
        double[][] motif;
        if( species >= 0 && model instanceof AbstractPhyloModel) {
            motif = ((AbstractPhyloModel) model).getCondProbsPerSpecies(species);
        } else {
            motif = model.getCondProbs();
        }
        
        double sumH = 0;
        for (int m = 0; m < motif.length; m++) {
			sumH += EvaluationUtil.determineIC(motif[m]);
        }
        System.out.println("Entropy:\t" + sumH);
    }
}
