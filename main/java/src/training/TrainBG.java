package training;

import java.util.Properties;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.models.AbstractModel;
import io.FileUtil;
import io.SampleUtil;
import models.AbstractAlignmentBasedModel;
import models.AbstractPhyloModel;
import models.ModelUtil;
import util.Config;

/**
 * This class manages the training of a {@link AbstractModel} (Background). It can handle {@link AbstractAlignmentBasedModel}s and
 * {@link AbstractPhyloModel}s.
 * 
 * @author Martin Nettling
 * 
 */
public class TrainBG {

	@SuppressWarnings("javadoc")
	public static void main(String[] args) throws Exception {
		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		String modelFileBG = Config.getProperty(props, "output.model.background", false).asString();

		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		// get datasets
		Sample[] samples = SampleUtil.getDataSets(props, alphabet);
		Sample sampleTrainBG = samples[2];

		// initialize new BG or load BG if possible (covered by #getNewBackground)
		AbstractModel newBG = ModelUtil.getNewBackground(props, alphabet);

		System.out.println("Training BG");
		System.out.println("Train BG on " + sampleTrainBG.getNumberOfElements() + " sequences");
		ModelUtil.trainModel(newBG, sampleTrainBG, 1e-3, 0, null); // train BG-Model on background

		FileUtil.writeFile(modelFileBG, newBG.toXML().toString());
	}
}
