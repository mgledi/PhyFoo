package training;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;
import java.util.Random;

import com.google.common.base.Joiner;
import com.google.common.base.Strings;

import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.termination.IterationCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.RecursiveSequence;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import io.FileUtil;
import io.PhyloSample;
import io.SampleUtil;
import models.AbstractAlignmentBasedModel;
import models.AbstractPhyloModel;
import models.AlignmentBasedModel;
import models.ModelUtil;
import models.PhyloBackground;
import models.PhyloBayesModel;
import models.PhyloPreparedAbstractModel;
import projects.dispom.PFMComparator;
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

		String TF = Config.getProperty(props, "TF", false).toString();
		String pValTable = Config.getProperty(props, "output.pValTable", "").toString();
		PhyloBayesModel.ASSUME_FILTERED_DATA = Config.getProperty(props, "model.fg.learnTemperature", "0").asInt();
		PhyloBackground.ASSUME_FILTERED_DATA = Config.getProperty(props, "model.bg.learnTemperature", "0").asInt();
		PhyloPreparedAbstractModel.MOTIF_QUALITY_THRESHOLD = Config.getProperty(props, "algorithm.training.motif_qualitiy_threshold", "0.0").asDouble();
		PhyloPreparedAbstractModel.PROB_THRESH_FOR_WEIGHTS = Config.getProperty(props, "algorithm.training.prob_thresh_for_weights", "1.1").asDouble();

		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		// get trained SHM to predict positions
		SingleHiddenMotifMixture bestSHM = ModelUtil.getTrainedSHM(props);
		sbR.append("TF[['" + TF + "']][['all']] = ");
		sbR.append(R.linearizeMatrix(PWMUtil.getCondProbs(bestSHM.model[0], false)) + ";\n");
		// get datasets
		PhyloSample sampleTrainFG = SampleUtil.getDataSetFGTrain(props, alphabet);

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
		AlignmentBasedModel motif = new AlignmentBasedModel(length, ((StrandModel) bestSHM.model[0]).model[0].getMaximalMarkovOrder(),
		        bestSHM.getAlphabetContainer());
		// fill species specific sequence containers and corresponding weights
		Sequence<?>[][] speciesSeqs = new Sequence<?>[species.size()][];
		double[][] speciesWeights = new double[species.size()][];
		for (int speciesIdx = 0; speciesIdx < species.size(); speciesIdx++) {
			parentToContent.clear();
			speciesWeights[speciesIdx] = new double[weights.length];
			speciesSeqs[speciesIdx] = new Sequence[weights.length];
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
					speciesWeights[speciesIdx][k] = weights[i];
					speciesSeqs[speciesIdx][k] = speciesSeq;
					k++;
				} else {
					ignoredSequences[speciesIdx]++;
				}
			}
			// copy arrays to avoid null values
			speciesWeights[speciesIdx] = Arrays.copyOf(speciesWeights[speciesIdx], k);
			speciesSeqs[speciesIdx] = Arrays.copyOf(speciesSeqs[speciesIdx], k);
		}
		StringBuilder sequenceStats = new StringBuilder();
		sequenceStats.append("Used sequences:\n");
		sequenceStats.append(Joiner.on(" & ").join(species)).append(" \\\\ \\hline \n");
		sequenceStats.append(TF + " & ");
		for (int i = 0; i < ignoredSequences.length; i++) {
			sequenceStats.append(speciesSeqs[0].length - ignoredSequences[i]).append(" & ");
		}
		sequenceStats.append(" \\\\ \\hline \n");
		sequenceStats.append(TF + " & ");
		for (int i = 0; i < ignoredSequences.length; i++) {
			sequenceStats.append((speciesSeqs[0].length - ignoredSequences[i]) * 100. / speciesSeqs[0].length).append(" & ");
		}
		sequenceStats.append(" \\\\ \\hline \n");
		sequenceStats.append(TF + " & ");
		for (int i = 0; i < ignoredSequences.length; i++) {
			sequenceStats.append(Util.round(Util.sum(speciesWeights[i]), 4)).append(" & ");
		}
		sequenceStats.append(" \\\\ \\hline \n");
		sequenceStats.append(TF + " & ");
		double refSum = Util.sum(speciesWeights[0]);
		for (int i = 0; i < ignoredSequences.length; i++) {
			sequenceStats.append(Util.round(Util.sum(speciesWeights[i]) / refSum * 100, 4)).append(" & ");
		}
		sequenceStats.append(" \\\\ \\hline \n");
		System.out.println(sequenceStats);

		// calculate and output spescies specific motifs
		System.out.println("Calculate species-specific motifs");
		for (int speciesIdx = 0; speciesIdx < species.size(); speciesIdx++) {
			// initialize new SHM
			AlignmentBasedModel.SKIP_TRAINING_STEPS = 0;
			motif.train(new Sample("Species sample", speciesSeqs[speciesIdx]), speciesWeights[speciesIdx]);
			sbR.append("TF[['" + TF + "']][['" + species.get(speciesIdx) + "']] = ");
			sbR.append(R.linearizeMatrix(PWMUtil.getCondProbs(motif, false)) + ";\n");
		}
		System.out.println(sbR);

		if (!pValTable.isEmpty()) {
			System.out.println("Determine P-Values");
			if (Config.getProperty(props, "species1.index", "-1").asInt() > -1 && Config.getProperty(props, "species2.index", "-1").asInt() > -1) {
				int s1 = Config.getProperty(props, "species1.index", false).asInt(), s2 = Config.getProperty(props, "species2.index", false).asInt();
				double val = Util.round(calculatePValue(speciesSeqs, speciesWeights, s1, s2, motif), 4);
				FileUtil.writeFile(pValTable + "_" + s1 + "_" + s2, val + "");
			} else {
				double[][] vals = new double[species.size()][species.size()];
				for (int s1 = 0; s1 < species.size(); s1++) {
					Arrays.fill(vals[s1], -1);
					for (int s2 = s1 + 1; s2 < species.size(); s2++) {
						System.out.println("Calculate pVal for motifs of species " + species.get(s1) + " and species " + species.get(s2));
						vals[s1][s2] = Util.round(calculatePValue(speciesSeqs, speciesWeights, s1, s2, motif), 4);
					}
				}
				String table = toLatex(vals, species, TF, motif.getMaximalMarkovOrder());
				System.out.println(table);
				FileUtil.writeFile(pValTable, sequenceStats + "\n\n" + table);
			}
		} else {
			System.out.println("Skip calculating pvalues.");
		}
	}

	/**
	 * Generates a latex table from the given vals and species.
	 */
	static String toLatex(double[][] vals, List<String> species, String TF, byte order) {
		StringBuilder body = new StringBuilder();
		body.append(" & ").append(Joiner.on(" & ").join(species)).append(" \\\\ \\hline \n");
		for (int i = 0; i < vals.length; i++) {
			body.append(species.get(i) + " & ");
			for (int j = 0; j < vals.length; j++) {
				body.append(vals[i][j] == -1 ? "" : vals[i][j] < 0.05 ? "\\textbf{" + vals[i][j] + "}" : vals[i][j])
				        .append(j < vals.length - 1 ? " & " : " \\\\ \\hline \n");
			}
		}
		StringBuilder table = new StringBuilder();
		table.append("\\begin{table}[!ht] \n");
		table.append("\\center");
		table.append("{\\footnotesize \n");
		table.append("\\begin{tabular}{|l").append(Strings.repeat("|c", species.size())).append("|} \\hline \n");
		table.append(body);
		table.append("\\end{tabular} \n");
		table.append("} \n");
		table.append("\\caption{p-Values for hypothesis that two species-specific motifs are generated by the same inhomogeous markov model of order "
		        + order + " for " + TF + "}\n");
		table.append("\\end{table} \n");

		return table.toString();
	}

	static double calculatePValue(Sequence<?>[][] speciesSeqs, double[][] speciesWeights, int s1, int s2, AlignmentBasedModel modelContainer)
	        throws EmptySampleException, WrongAlphabetException, Exception {
		PFMComparator.SymmetricKullbackLeiblerDivergence klMeasure = new PFMComparator.SymmetricKullbackLeiblerDivergence(0);
		double[][] pfm1, pfm2;
		modelContainer.train(new Sample("Species sample", speciesSeqs[s1]), speciesWeights[s1]);
		pfm1 = Util.arraycopy(modelContainer.getCondProbs());
		modelContainer.train(new Sample("Species sample", speciesSeqs[s2]), speciesWeights[s2]);
		pfm2 = Util.arraycopy(modelContainer.getCondProbs());
		double klOrig = klMeasure.getDistance(pfm1, pfm2, 0);

		// start test
		int countSmaller = 0, tests = 1000;
		for (int i = 0; i < tests; i++) {
			Sequence<?>[] speciesSeqs1 = Arrays.copyOf(speciesSeqs[s1], speciesSeqs[s1].length);
			double[] speciesWeights1 = Arrays.copyOf(speciesWeights[s1], speciesWeights[s1].length);
			Sequence<?>[] speciesSeqs2 = Arrays.copyOf(speciesSeqs[s2], speciesSeqs[s2].length);
			double[] speciesWeights2 = Arrays.copyOf(speciesWeights[s2], speciesWeights[s2].length);

			labelSwitching(speciesSeqs1, speciesWeights1, speciesSeqs2, speciesWeights2);
			modelContainer.train(new Sample("Species sample 1", speciesSeqs1), speciesWeights1);
			pfm1 = Util.arraycopy(modelContainer.getCondProbs());
			modelContainer.train(new Sample("Species sample 2", speciesSeqs2), speciesWeights2);
			pfm2 = Util.arraycopy(modelContainer.getCondProbs());
			double klTest = klMeasure.getDistance(pfm1, pfm2, 0);

			if (klTest <= klOrig) {
				countSmaller++;
			}
		}
		return countSmaller * 1. / tests;
	}

	static void labelSwitching(Sequence<?>[] speciesSeqs1, double[] speciesWeights1, Sequence<?>[] speciesSeqs2, double[] speciesWeights2) {
		Random rand = new Random();

		int offset = speciesSeqs1.length;
		for (int i = 0; i < speciesSeqs1.length; i++) {
			int newIdx = rand.nextInt(speciesSeqs1.length + speciesSeqs2.length);
			Sequence<?> tmpSeq = speciesSeqs1[i];
			double tmpWeight = speciesWeights1[i];
			// no label switching
			if (newIdx < speciesSeqs1.length) {
				speciesSeqs1[i] = speciesSeqs1[newIdx];
				speciesSeqs1[newIdx] = tmpSeq;
				speciesWeights1[i] = speciesWeights1[newIdx];
				speciesWeights1[newIdx] = tmpWeight;
			} else {
				speciesSeqs1[i] = speciesSeqs2[newIdx - offset];
				speciesSeqs2[newIdx - offset] = tmpSeq;
				speciesWeights1[i] = speciesWeights2[newIdx - offset];
				speciesWeights2[newIdx - offset] = tmpWeight;
			}
		}
		for (int i = 0; i < speciesSeqs2.length; i++) {
			int newIdx = rand.nextInt(speciesSeqs1.length + speciesSeqs2.length);
			Sequence<?> tmpSeq = speciesSeqs2[i];
			double tmpWeight = speciesWeights2[i];
			// no label switching
			if (newIdx >= speciesSeqs1.length) {
				speciesSeqs2[i] = speciesSeqs2[newIdx - offset];
				speciesSeqs2[newIdx - offset] = tmpSeq;
				speciesWeights2[i] = speciesWeights2[newIdx - offset];
				speciesWeights2[newIdx - offset] = tmpWeight;
			} else {
				speciesSeqs2[i] = speciesSeqs1[newIdx];
				speciesSeqs1[newIdx] = tmpSeq;
				speciesWeights2[i] = speciesWeights1[newIdx];
				speciesWeights1[newIdx] = tmpWeight;
			}
		}
	}
}
