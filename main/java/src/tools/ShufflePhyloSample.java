package tools;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Properties;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import io.PhyloSample;
import io.SampleUtil;
import util.Config;

public class ShufflePhyloSample {
	public static void main(String[] args) throws IllegalArgumentException, EmptySampleException, WrongAlphabetException, WrongLengthException, IOException,
	        WrongSequenceTypeException, DoubleSymbolException {
		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		String inputFile = Config.getProperty(props, "input.file", false).asString();
		String outputFile = Config.getProperty(props, "output.file", false).asString();
		String newick = Config.getProperty(props, "model.newick", false).asString();

		PhyloSample sample = SampleUtil.getDataSet(inputFile, newick, new AlphabetContainer(new UnobservableDNAAlphabet()));
		PhyloSample shuffledSample = sample.shuffle(newick);

		String fasta = SampleUtil.toFasta(shuffledSample);
		File output = new File(outputFile);
		FileWriter fw = new FileWriter(output);
		fw.write(fasta);
		fw.close();
	}
}
