package training;

import java.util.Properties;

import org.apache.log4j.Logger;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.models.AbstractModel;
import io.FileUtil;
import io.PhyloSample;
import io.SampleUtil;
import models.AbstractAlignmentBasedModel;
import models.AbstractPhyloModel;
import models.ModelUtil;
import util.Config;

/**
 * This class manages the training of a {@link AbstractModel} (Background). It can handle {@link AbstractAlignmentBasedModel}s and {@link AbstractPhyloModel}s.
 * The resulting model is written to a defined XML-file.
 * 
 * @author Martin Nettling
 * 
 */
public class TrainBG {
	private static Logger LOGGER = Logger.getLogger(TrainBG.class);

	@SuppressWarnings("javadoc")
	public static void main(String[] args) throws Exception {
		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		String modelFileBG = Config.getProperty(props, "output.model.background", false).asString();

		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		// get datasets
		PhyloSample sampleTrainBG = SampleUtil.getDataSetBGTrain(props, alphabet);
		// initialize new BG or load BG if possible (covered by #getNewBackground)
		AbstractModel newBG = ModelUtil.getNewBackground(props, alphabet);

		LOGGER.info("Training BG");
		LOGGER.info("Train BG on " + sampleTrainBG.getNumberOfElements() + " alignments.");
		TrainingUtil.trainModel(newBG, sampleTrainBG, 1e-3, 0, null); // train BG-Model on background
		LOGGER.info("Writing model to " + modelFileBG);
		FileUtil.writeFile(modelFileBG, newBG.toXML().toString());
	}
}
