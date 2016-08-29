package training;

import org.apache.log4j.Logger;

import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.algorithms.optimization.termination.TerminationCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.models.Model;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import models.PhyloBackground;
import util.Util;

public class TrainingUtil {
	private static Logger LOGGER = Logger.getLogger(TrainingUtil.class);

	public static void trainModel(Model model, Sample data, double eps, int minSteps, double[] weights) throws Exception {
		byte order = -1;
		if (model instanceof StrandModel) {
			order = ((de.jstacs.models.AbstractModel) ((StrandModel) model).model[0]).getMaximalMarkovOrder();
		} else {
			order = ((de.jstacs.models.AbstractModel) model).getMaximalMarkovOrder();
			LOGGER.info("Start training " + model.getInstanceName() + "(" + order + ")");
		}
		double lastLog = Integer.MIN_VALUE, curLog = Integer.MIN_VALUE;

		int i = 0;
		while (i == 0 || i < minSteps || Math.abs(lastLog - curLog) > eps && i < 20) {
			if (weights == null) {
				model.train(data);
			} else {
				model.train(data, weights);
			}
			lastLog = curLog;
			curLog = Util.sum(model.getLogProbFor(data));
			LOGGER.info("Step: " + i + " | Difference = " + (curLog - lastLog));
			i++;
		}
	}

	/**
	 * 
	 * @param trainingFile
	 * @param shm
	 * @throws Exception
	 */
	public static void trainZoopsAndStrandParameters(String trainingFile, SingleHiddenMotifMixture shm) throws Exception {
		if (!((shm.model[0] instanceof StrandModel) && ((StrandModel) shm.model[0]).model[0] instanceof BayesianNetworkModel
		        || (shm.model[0] instanceof BayesianNetworkModel))) {
			throw new Exception("Expected BayesianNetworkModel as Motifmodel in the given SingleHiddenMotifMixture");
		}
		// remembering old TerminationCondition
		TerminationCondition tcTmp = shm.tc;
		// Overwriting TerminationCondition to do only one Iteration
		shm.tc = new IterationCondition(0);
		BayesianNetworkModel.SKIP_TRAINING = true;
		AlphabetContainer con = shm.getAlphabetContainer();
		SplitSequenceAnnotationParser parser = new SplitSequenceAnnotationParser("=", ";");
		Sample fgTrain = new Sample(con, new SparseStringExtractor(trainingFile, '>', parser));
		shm.setShiftCorrection(false);
		shm.train(fgTrain);
		BayesianNetworkModel.SKIP_TRAINING = false;

		// resetting TerminationCondition
		shm.tc = tcTmp;
	}

	public static String learnTree(Sample data, String initialTree) throws Exception {
		LOGGER.info("Learning Tree, starting at " + initialTree);
		PhyloBackground bg = new PhyloBackground((byte) 0, initialTree, data.getAlphabetContainer());
		PhyloBackground.EDGE_LEARNING = true;
		trainModel(bg, data, 1e-1, 0, null);
		PhyloBackground.EDGE_LEARNING = false;
		PhyloBackground.EDGE_LEARNING_ONE_LENGTH = false;
		return bg.getBNH().getVirtualTree(0).getNewickString();
	}
}
