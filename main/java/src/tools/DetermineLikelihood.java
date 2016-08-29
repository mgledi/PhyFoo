package tools;

import java.util.Properties;

import de.jstacs.data.Sample;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import io.SampleUtil;
import models.ModelUtil;
import models.PhyloBackground;
import models.PhyloBayesModel;
import models.PhyloPreparedAbstractModel;
import util.Config;
import util.Util;

public class DetermineLikelihood {

    public static void main(String[] args) throws Exception {
        Properties props = new Properties();
        Config.parseProperties(props, args, false);
        Config.replaceInternalPlaceholder(props);

        // get trained SHM to predict positions
        SingleHiddenMotifMixture trainedSHM = null;
        if (props.containsKey("input.model.shm") && props.getProperty("input.model.shm") != null && !props.getProperty("input.model.shm").isEmpty()) {
            System.out.println("====== Load Learned SHM from " + props.getProperty("input.model.shm") + " =====");
            trainedSHM = ModelUtil.getTrainedSHM(props);
        }
        PhyloPreparedAbstractModel trainedBackground = null;
        if (props.containsKey("input.model.background") && props.getProperty("input.model.background") != null && !props.getProperty("input.model.background").isEmpty()) {
            System.out.println("====== Load Learned Backgroundfrom " + props.getProperty("input.model.background") + " =====");
            trainedBackground = ModelUtil.getTrainedBackground(props);
        }
        // get datasets
        Sample[] samples = SampleUtil.getDataSets(props, trainedSHM.alphabets);

        PhyloBayesModel.ENABLE_CACHING = false;
        PhyloBackground.ENABLE_CACHING = false;

        // ##########################################################################
        // load learned model

        Sample sampleTrainFG = samples[0];
        Sample sampleTestFG = samples[1];
        Sample sampleTrainBG = samples[2];
        Sample sampleTestBG = samples[3];

        if (trainedSHM != null) {
            System.out.println("LogLikelihood FG FG Train:\t" + Util.sum(trainedSHM.getLogProbFor(sampleTrainFG)));
            System.out.println("LogLikelihood FG FG Test:\t" + Util.sum(trainedSHM.getLogProbFor(sampleTestFG)));

            System.out.println("LogLikelihood FG BG Train:\t" + Util.sum(trainedSHM.getLogProbFor(sampleTrainBG)));
            System.out.println("LogLikelihood FG BG Test:\t" + Util.sum(trainedSHM.getLogProbFor(sampleTestBG)));
        }

        if (trainedBackground != null) {
            System.out.println("LogLikelihood BG FG Train:\t" + Util.sum(trainedBackground.getLogProbFor(sampleTrainFG)));
            System.out.println("LogLikelihood BG FG Test:\t" + Util.sum(trainedBackground.getLogProbFor(sampleTestFG)));

            System.out.println("LogLikelihood BG BG Train:\t" + Util.sum(trainedBackground.getLogProbFor(sampleTrainBG)));
            System.out.println("LogLikelihood BG BG Test:\t" + Util.sum(trainedBackground.getLogProbFor(sampleTestBG)));
        }
    }
}
