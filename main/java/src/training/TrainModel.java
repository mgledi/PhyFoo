package training;

import java.util.ArrayList;
import java.util.Properties;
import java.util.concurrent.TimeUnit;

import org.apache.log4j.Logger;

import algorithm.EMHandler;
import de.jstacs.algorithms.optimization.termination.CombinedCondition;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.algorithms.optimization.termination.TimeCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import io.FileUtil;
import io.PhyloSample;
import io.SampleUtil;
import models.AbstractAlignmentBasedModel;
import models.AbstractPhyloModel;
import models.ModelUtil;
import models.PhyloPreparedAbstractModel;
import util.Config;
import util.PWMUtil;
import util.Util;

/**
 * This class manages the training of a {@link SingleHiddenMotifMixture}. The {@link SingleHiddenMotifMixture} can encapsulate
 * {@link AbstractAlignmentBasedModel}s and {@link AbstractPhyloModel}s. The resulting model is written to a defined XML-file.
 * 
 * @author Martin Nettling
 */
public class TrainModel {
	private static Logger LOGGER = Logger.getLogger(TrainModel.class);

	@SuppressWarnings("javadoc")
	public static void main(String[] args) throws Exception {
		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		String modelFileSHM = Config.getProperty(props, "output.model.shm", false).asString();

		if (Config.getProperty(props, "algorithm.continue", "false").asBoolean()) {
			props.setProperty("input.model.shm", modelFileSHM);
		}

		PhyloPreparedAbstractModel.MOTIF_QUALITY_THRESHOLD = Config.getProperty(props, "algorithm.training.motif_qualitiy_threshold", "0.0").asDouble();
		PhyloPreparedAbstractModel.PROB_THRESH_FOR_WEIGHTS = Config.getProperty(props, "algorithm.training.prob_thresh_for_weights", "0.8").asDouble();
		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		PhyloSample sampleTrainFG = SampleUtil.getDataSetFGTrain(props, alphabet);

		CombinedCondition stopCondition = new CombinedCondition(3, new IterationCondition(Config.getProperty(props, "algorithm.emsteps", "100").asInt()),
		        new SmallDifferenceOfFunctionEvaluationsCondition(Config.getProperty(props, "algorithm.tc.smalldifference", "0.0001").asDouble()),
		        new TimeCondition(TimeUnit.MINUTES.toSeconds(Config.getProperty(props, "algorithm.maxrunningtime", "180").asInt())));

		// initialize new SHM or load SHM if possible (covered by #getNewSHM)
		SingleHiddenMotifMixture newSHM = ModelUtil.getNewSHM(props, alphabet, stopCondition);

		LOGGER.info("Training BG");
		LOGGER.info("Train FG on " + sampleTrainFG.getNumberOfElements() + " sequences");
		// Train flanking model if it is not trained yet
		if (Config.getProperty(props, "algorithm.continue", "false").asBoolean()) {
			LOGGER.info("Skip training of flanking model");
		} else {
			TrainingUtil.trainModel(newSHM.model[1], sampleTrainFG, 1e-3, 0, null); // train flanking Model on foreground
		}

		// TODO: Allow configuring Tree Learning
		if (Config.getProperty(props, "algorithm.pretrainonknownpositions", "false").asBoolean()) {
			// generate sample, containing only the true motifs contained in dataFG
			int[] truePositionsSet = SampleUtil.getTruePositionSet(sampleTrainFG, "MOTIF_POS");
			ArrayList<Sequence<?>> motifs = new ArrayList<Sequence<?>>(truePositionsSet.length);
			for (int i = 0; i < truePositionsSet.length; i++) {
				// if no motif was found, don't add it to the list of training motifs
				if (truePositionsSet[i] != -1) {
					motifs.add(sampleTrainFG.getElementAt(i).getSubSequence(truePositionsSet[i], newSHM.model[0].getLength()));
				}
			}
			LOGGER.info("Start PreTraining on " + motifs.size() + " motifs.");
			Sample motifSample = new Sample("True Motifs", motifs.toArray(new Sequence[0]));
			TrainingUtil.trainModel(newSHM.model[0], motifSample, 1e-4, 2, null);
			LOGGER.info("Finished Training. Motif learned:");
			Util.visualizeMatrix(PWMUtil.getCondProbs(newSHM.model[0], true));
			LOGGER.info("===============");
		}

		newSHM.setOutputStream(System.out);

		// fill EMHandler
		EMHandler emhandler = EMHandler.train(newSHM, Config.getProperty(props, "algorithm.emrestarts", "1").asInt(), sampleTrainFG, modelFileSHM);

		SingleHiddenMotifMixture bestModel = emhandler.getBestModel();

		LOGGER.info("Best Model: " + bestModel.getScoreForBestRun());
		LOGGER.info("Writing model to " + modelFileSHM);
		FileUtil.writeFile(modelFileSHM, bestModel.toXML().toString());
	}
}

