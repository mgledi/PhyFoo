package tools;

import java.util.Properties;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.motifDiscovery.MotifDiscoverer.KindOfProfile;
import io.SampleUtil;
import models.ModelUtil;
import util.Config;

/**
 * --basedir="e:/Dropbox/promotion/edgeLengthTest_FS81alpha/biological_tests/CTCF_gapless" <br>
 * --input.dataset.fg=%basedir%/dataset/data_1.0.fg <br>
 * --input.dataset.bg=%basedir%/dataset/data_1.0.bg <br>
 * --model.newick=(SPECIES_0:1.0,SPECIES_1:1.0,SPECIES_2:1.0,SPECIES_3:1.0,SPECIES_4:1.0) <br>
 * --input.splitseed=0 <br>
 * --output.dataset.fg.train=%basedir%/dataset/data_1.0.fg.train <br>
 * --output.dataset.fg.test=%basedir%/dataset/data_1.0.fg.test <br>
 * --output.dataset.bg.train=%basedir%/dataset/data_1.0.bg.train <br>
 * --output.dataset.bg.test=%basedir%/dataset/data_1.0.bg.test <br>
 * </code>
 */
public class ExtractBindingSitesPerSpecies {
	public static void main(String[] args) throws Exception {
		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		String species = props.getProperty("output.species");

		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		// get datasets
		Sample[] samples = SampleUtil.getDataSets(props, alphabet);
		Sample sampleTrainFG = samples[0], sampleTestFG = samples[1], sampleTrainBG = samples[2], sampleTestBG = samples[3];

		SingleHiddenMotifMixture shm = ModelUtil.getTrainedSHM(props);
//		SingleHiddenMotifMixture shm = ModelUtil.getNewSHM(props, alphabet, new IterationCondition(1));
//		shm.train(sampleTrainFG);
		
		// create the corresponding motifSample, this sample is needed, when only maximization is done
		Sample sampleFG = Sample.union(sampleTrainFG, sampleTestFG);
		for (int i = 0; i < sampleFG.getNumberOfElements(); i++) {
			// calculate probabilities for of shm for each position in the current sequence
			// here seems to be bug, we must call getProfileOfScores first
			shm.getProfileOfScoresFor(0, 0, sampleFG.getElementAt(i), 0, KindOfProfile.NORMALIZED_CONDITIONAL);
			
			// extract binding site
			Sequence<?> bindingSite = sampleFG.getElementAt(i).getSubSequence(0,shm.model[0].getLength());
			// determine most probable strand
			if(shm.model[0] instanceof StrandModel) {
				double p0 = ((StrandModel)shm.model[0]).getLogProbFor(0, bindingSite);
				double p1 = ((StrandModel)shm.model[0]).getLogProbFor(1, bindingSite);
				
				MultiDimensionalDiscreteSequence parent = (MultiDimensionalDiscreteSequence) sampleFG.getElementAt(i);
				for(int o=0; o < parent.getNumberOfSequences(); o++) {
					for(SequenceAnnotation sa : parent.getSequence(o).getAnnotation()) {
						if(sa.getType().equals("species")) {
							species = sa.getIdentifier();
						}
					}
					System.out.println(parent.getSequence(o).getSubSequence(0,shm.model[0].getLength()) + "\t" + p0 + "\t" + 0 + "\t" + species);
					System.out.println(parent.getSequence(o).getSubSequence(0,shm.model[0].getLength()).reverseComplement() + "\t" + p1 + "\t" + 1 + "\t" + species);
				}
				System.out.println("---------");
				
			}
		}
	}
}
