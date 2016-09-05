package training;

import java.util.Properties;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import classification.ClassificationUtil;
import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TimeCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.models.AbstractModel;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import io.FileUtil;
import io.SampleUtil;
import models.AbstractAlignmentBasedModel;
import models.AbstractPhyloModel;
import models.ModelUtil;
import models.PhyloBackground;
import models.PhyloBayesModel;
import models.PhyloPreparedAbstractModel;
import util.Config;
import util.PWMUtil;
import util.Util;

/**
 * This class contains methods to refine an already trained {@link SingleHiddenMotifMixture}. It needs a complete Model
 * configuration and a xml file containing the model to be refined. The {@link SingleHiddenMotifMixture} is allowed to
 * consist of {@link AbstractAlignmentBasedModel} and {@link AbstractPhyloModel}.
 * 
 * @author Martin Nettling
 * 
 */
public class PostTrainModel {

	private static Logger LOGGER = Logger.getLogger(PostTrainModel.class);

    @SuppressWarnings("javadoc")
    public static void main(String[] args) throws Exception {
        // do classification
        PhyloBayesModel.ENABLE_CACHING = false;
        PhyloBayesModel.CLONING_ALLOWED = false;
        
        
        Properties props = new Properties();
        Config.parseProperties(props, args, false);
        Config.replaceInternalPlaceholder(props);

        PhyloBayesModel.ASSUME_FILTERED_DATA = Config.getProperty(props, "model.fg.learnTemperature", "0").asInt();
        PhyloBackground.ASSUME_FILTERED_DATA = Config.getProperty(props, "model.bg.learnTemperature", "0").asInt(); 
    	PhyloPreparedAbstractModel.MOTIF_QUALITY_THRESHOLD = Config.getProperty(props, "algorithm.training.motif_qualitiy_threshold","0.0").asDouble();
    	PhyloPreparedAbstractModel.PROB_THRESH_FOR_WEIGHTS = Config.getProperty(props, "algorithm.training.prob_thresh_for_weights","0.8").asDouble();
    	
        AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

        // read needed variables to find wrong configuration early
        String resultFile = Config.getProperty(props, "output.results", false).asString();
        String modelFileSHM = Config.getProperty(props, "output.model.shm", false).asString();
        String modelFileBG = Config.getProperty(props, "output.model.background", false).asString();

        // get trained SHM to predict positions
        SingleHiddenMotifMixture bestSHM = ModelUtil.getTrainedSHM(props);
        PhyloPreparedAbstractModel bestBG = ModelUtil.getTrainedBackground(props);

        // get datasets
        Sample[] samples = SampleUtil.getDataSets(props, bestSHM.alphabets);
        Sample sampleTrainFG = samples[0], sampleTestFG = samples[1], sampleTrainBG = samples[2], sampleTestBG = samples[3];
        // get predictions on training foreground dataset
        double weights[] = bestSHM.getSeqWeigths(sampleTrainFG);

        // create the corresponding motifSample, this sample is needed, when only maximization is done
        Sample motifSample = SampleUtil.generateTrainingSample(sampleTrainFG, bestSHM.model[0].getLength());

        // initialize new SHM
        CombinedCondition stopCondition = new CombinedCondition(3,
                new IterationCondition(Config.getProperty(props, "algorithm.emsteps", "100").asInt()),
                new SmallDifferenceOfFunctionEvaluationsCondition(Config.getProperty(props, "algorithm.tc.smalldifference","0.0001").asDouble()),
                new TimeCondition(TimeUnit.MINUTES.toSeconds(Config.getProperty(props, "algorithm.maxrunningtime","180").asInt())));

        SingleHiddenMotifMixture newSHM = ModelUtil.getNewSHM(props, alphabet, stopCondition);
        AbstractModel trainingBG = ModelUtil.getNewBackground(props, alphabet);

        // only SHM and StrandModel are cloned. references to not clone supporting models stay
        SingleHiddenMotifMixture trainingSHM = (SingleHiddenMotifMixture) bestSHM.clone();
		trainingSHM.tc = stopCondition;
        // use empty phylo models from initial SHM. (new references, bestSHM is untouched)
        if(trainingSHM.model[0] instanceof StrandModel) {
			((StrandModel) trainingSHM.model[0]).model[0] = ((StrandModel) newSHM.model[0]).model[0];
        }
        trainingSHM.model[1] = newSHM.model[1];

        // copy conditional probabilities to training SHM
		double[][] condProbs = PWMUtil.getCondProbs(bestSHM.model[0], false);
        if(trainingSHM.model[0] instanceof StrandModel) {
            ((PhyloPreparedAbstractModel) ((StrandModel)trainingSHM.model[0]).model[0]).setCondProbs(condProbs);
        } else {
            ((PhyloPreparedAbstractModel) trainingSHM.model[0]).setCondProbs(condProbs);
        }
        
        // copy probabilities from old bg to new bg models
        ((PhyloPreparedAbstractModel)trainingBG).setCondProbs(bestBG.getCondProbs());
        ((PhyloPreparedAbstractModel)trainingSHM.model[1]).setCondProbs(((PhyloPreparedAbstractModel)bestSHM.model[1]).getCondProbs());
        
        // train flanking model and background model
		TrainingUtil.trainModel(trainingBG, sampleTrainBG, 1e-3, 0, null); // train BG-Model on background
		TrainingUtil.trainModel(trainingSHM.model[1], sampleTrainFG, 1e-3, 0, null); // train flanking Model on foreground
        
        // train on extracted gammas
        if (Config.getProperty(props, "algorithm.maximization_only", "true").asBoolean()) {
			TrainingUtil.trainModel(trainingSHM.model[0], motifSample, 1e-4, 2, weights);
        } else { 
            trainingSHM.train(sampleTrainFG);
        }
        
		LOGGER.info("Writing SHM model to " + modelFileSHM);
		FileUtil.writeFile(modelFileSHM, trainingSHM.toXML().toString());
		LOGGER.info("Writing BG model to " + modelFileBG);
		FileUtil.writeFile(modelFileBG, trainingBG.toXML().toString());
        // perform classification tests
        

        PhyloBayesModel.ENABLE_CACHING = false;
        PhyloBayesModel.CLONING_ALLOWED = false;
		LOGGER.info("Start tests on " + sampleTestFG.getNumberOfElements() + " FG-seqs and " + sampleTestBG.getNumberOfElements() + " BG-seqs.");
        double[] classTest = ClassificationUtil
                .performClassificationTest(newSHM, trainingBG, sampleTestFG, sampleTestBG);

        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < classTest.length; i++) {
            sb.append(Util.round(classTest[i], 4) + "\n");
        }
        
        // store model and results
		LOGGER.info("Writing results to " + resultFile);
		FileUtil.writeFile(resultFile, sb.toString());

    }
}
