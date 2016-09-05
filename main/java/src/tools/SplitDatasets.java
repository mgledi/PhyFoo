package tools;

import java.io.IOException;
import java.util.Properties;

import org.apache.log4j.Logger;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.Sample;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import io.FileUtil;
import io.SampleUtil;
import util.Config;

/**
 * This tool splits an input dataset according to given seed into foreground-train, foreground-test, background-train,
 * and background-test.
 * 
 * Example Parameterization<br>
 * <code>
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
 * 
 * @author Martin Nettling
 * 
 */
public class SplitDatasets {
	private static Logger LOGGER = Logger.getLogger(SplitDatasets.class);
    @SuppressWarnings("javadoc")
    public static void main(String[] args)
            throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException,
            WrongSequenceTypeException, IllegalArgumentException, DoubleSymbolException {
        Properties props = new Properties();
        Config.parseProperties(props, args, false);
        Config.replaceInternalPlaceholder(props);

        AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

        // get datasets
        Sample[] samples = SampleUtil.getDataSets(props, alphabet);

		LOGGER.info("Writing fasta for positive training dataset to " + Config.getProperty(props, "output.dataset.fg.train", false).toString());
		FileUtil.writeFile(Config.getProperty(props, "output.dataset.fg.train", false).toString(), SampleUtil.toFasta(samples[0]));
		LOGGER.info("Writing fasta for positive testing dataset to " + Config.getProperty(props, "output.dataset.fg.test", false).toString());
		FileUtil.writeFile(Config.getProperty(props, "output.dataset.fg.test", false).toString(), SampleUtil.toFasta(samples[1]));
		LOGGER.info("Writing fasta for negative traing dataset to " + Config.getProperty(props, "output.dataset.bg.train", false).toString());
		FileUtil.writeFile(Config.getProperty(props, "output.dataset.bg.train", false).toString(), SampleUtil.toFasta(samples[2]));
		LOGGER.info("Writing fasta for negative testing dataset to " + Config.getProperty(props, "output.dataset.bg.test", false).toString());
		FileUtil.writeFile(Config.getProperty(props, "output.dataset.bg.test", false).toString(), SampleUtil.toFasta(samples[3]));
    }
}
