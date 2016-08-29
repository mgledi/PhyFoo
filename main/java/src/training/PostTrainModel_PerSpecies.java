package training;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Properties;

import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.RecursiveSequence;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import io.SampleUtil;
import models.AbstractAlignmentBasedModel;
import models.AbstractPhyloModel;
import models.AlignmentBasedModel;
import models.ModelUtil;
import models.PhyloBackground;
import models.PhyloBayesModel;
import models.PhyloPreparedAbstractModel;
import util.Config;
import util.PWMUtil;
import util.R;
import util.Util;

/**
 * This class contains methods to refine an already trained {@link SingleHiddenMotifMixture} per species. It needs a complete Model configuration and a xml file
 * containing the model to be refined. The {@link SingleHiddenMotifMixture} is allowed to consist of {@link AbstractAlignmentBasedModel} and
 * {@link AbstractPhyloModel}.
 * 
 * @author Martin Nettling
 * 
 */
public class PostTrainModel_PerSpecies {

	@SuppressWarnings("javadoc")
	public static void main(String[] args) throws Exception {
		// do classification
		StringBuilder sbR = new StringBuilder();
		PhyloBayesModel.ENABLE_CACHING = false;
		PhyloBayesModel.CLONING_ALLOWED = false;

		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		String TF= Config.getProperty(props, "TF",false).toString();
		PhyloBayesModel.ASSUME_FILTERED_DATA = Config.getProperty(props, "model.fg.learnTemperature", "0").asInt();
		PhyloBackground.ASSUME_FILTERED_DATA = Config.getProperty(props, "model.bg.learnTemperature", "0").asInt();
		PhyloPreparedAbstractModel.MOTIF_QUALITY_THRESHOLD = Config.getProperty(props, "algorithm.training.motif_qualitiy_threshold", "0.0").asDouble();
		PhyloPreparedAbstractModel.PROB_THRESH_FOR_WEIGHTS = Config.getProperty(props, "algorithm.training.prob_thresh_for_weights", "1.1").asDouble();

		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		// get trained SHM to predict positions
		SingleHiddenMotifMixture bestSHM = ModelUtil.getTrainedSHM(props);
		sbR.append("TF[['"+TF+"']][['all']] = ");
		sbR.append(R.linearizeMatrix(PWMUtil.getCondProbs(bestSHM.model[0], false)) + ";\n");
		// get datasets
		Sample[] samples = SampleUtil.getDataSets(props, alphabet);
		Sample sampleTrainFG = samples[0];

		AlignmentBasedModel.SKIP_TRAINING_STEPS = 2;
		// SingleHiddenMotifMixture bestSHM = ModelUtil.getNewSHM(props, alphabet, new IterationCondition(1));
		bestSHM.tc = new IterationCondition(1);
		bestSHM.train(sampleTrainFG);

		double weights[] = ((AlignmentBasedModel) ((StrandModel) bestSHM.model[0]).model[0]).lastWeights;
		Sample motifSample = ((AlignmentBasedModel) ((StrandModel) bestSHM.model[0]).model[0]).lastSequences;

		// ############################################################################################################################################
		// ############################################################################################################################################
		// cache for all parent sequences. Needed due to massive abuse of clone in JStacs.
		HashMap<MultiDimensionalDiscreteSequence, MultiDimensionalDiscreteSequence> parentToContent = new HashMap<MultiDimensionalDiscreteSequence, MultiDimensionalDiscreteSequence>();
		String oldNewick = Config.getPropertyWithFallBack(props, false, "model.bg.newick", "model.fg.newick", "model.newick").asString();
		ArrayList<String> species = Util.getOrderedArrayListFromTree(oldNewick);
		int ignoredSequences[] = new int[species.size()];
		int length = motifSample.getElementAt(0).getLength();
		AlignmentBasedModel motif = new AlignmentBasedModel(
				length, 
		        ((StrandModel) bestSHM.model[0]).model[0].getMaximalMarkovOrder(),
		        bestSHM.getAlphabetContainer());
		for (int speciesIdx = 0; speciesIdx < species.size(); speciesIdx++) {
			parentToContent.clear();
			double[] speciesWeights = new double[weights.length];
			Sequence<?>[] speciesSeqs = new Sequence[weights.length];
			System.out.println("Transform sample for species " + species.get(speciesIdx));
			int k = 0;
			for (int i = 0; i < weights.length; i++) {
				MultiDimensionalDiscreteSequence parent = (MultiDimensionalDiscreteSequence) ((RecursiveSequence<?>) motifSample.getElementAt(i)).getParent();
				if (!parentToContent.containsKey(parent)) {
					parentToContent.put(parent, new MultiDimensionalDiscreteSequence(null, (ByteSequence) parent.getSequence(speciesIdx)));
				}
				int start = ((RecursiveSequence<?>) motifSample.getElementAt(i)).getIndex(0);
				Sequence<?> speciesSeq = parentToContent.get(parent).getSubSequence(start, length);

				if (!speciesSeq.toString().contains("-")) {
					speciesWeights[k] = weights[i];
					speciesSeqs[k] = speciesSeq;
					k++;
				} else {
					ignoredSequences[speciesIdx]++;
				}
			}
			speciesWeights = Arrays.copyOf(speciesWeights, k);
			speciesSeqs = Arrays.copyOf(speciesSeqs, k);
			// initialize new SHM
			AlignmentBasedModel.SKIP_TRAINING_STEPS = 0;
			motif.train(new Sample("Species sample", speciesSeqs), speciesWeights);

			sbR.append("TF[['" + TF + "']][['" + species.get(speciesIdx) + "']] = ");
			sbR.append(R.linearizeMatrix(PWMUtil.getCondProbs(motif, false)) + ";\n");
		}
		System.out.println(sbR);
		System.out.println(Arrays.toString(ignoredSequences));
	}
}
