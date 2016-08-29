package io;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Properties;
import java.util.Random;
import java.util.logging.Logger;

import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.data.sequences.SimpleDiscreteSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.data.sequences.annotation.SplitSequenceAnnotationParser;
import de.jstacs.io.SparseStringExtractor;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import util.Config;
import util.Util;

public class SampleUtil {
	static Logger LOGGER = Logger.getLogger(SampleUtil.class.getSimpleName());
    /**
     * This method generates a set of TFBSs from the given sample according to
     * {@link SingleHiddenMotifMixture#setTrainData}
     * 
     * @param sample
     * @param motifLength
     * @return
     */
    public static Sample generateTrainingSample(Sample sample, int motifLength) throws EmptySampleException,
            WrongAlphabetException {
        ArrayList<Sequence<?>> seqs = new ArrayList<Sequence<?>>();
        for (int i = 0; i < sample.getNumberOfElements(); i++) {
            for (int u = 0; u <= sample.getElementAt(i).getLength() - motifLength; u++) {
                seqs.add(sample.getElementAt(i).getSubSequence(u, motifLength));
            }
        }
        return new Sample("Motif sample of " + sample.getAnnotation(), seqs.toArray(new Sequence[0]));
    }

    /**
     * This method loads a {@link PhyloSample} from the given dataFile and returns it.
     * 
     * @param dataFile
     * @param newick
     * @param con
     * @return a PhyloSample represented by the given dataFile.
     */
	public static PhyloSample getDataSet(String dataFile, String newick, AlphabetContainer con)
            throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException,
            WrongSequenceTypeException {
        ArrayList<String> speciesOrder = Util.getOrderedArrayListFromTree(newick);
        SplitSequenceAnnotationParser parser = new SplitSequenceAnnotationParser("=", ";");
        Sample rawDataFGtrain = new Sample(con, new SparseStringExtractor(dataFile, '>', parser));
		PhyloSample sampleTrainFG = new PhyloSample(rawDataFGtrain, speciesOrder);
        return sampleTrainFG;
    }

    /**
     * This method reads the dataFile into a {@link PhyloSample} and tries to split the Sample into a training and a
     * testing sample according to the given split seed. The method returns an array of two samples. The first sample
     * represents the training part, the second sample is the testing part.
     * 
     * @param dataFile
     * @param newick
     * @param con
     * @param trainingSeqs
     * @param testingSeqs
     * @param splitseed
     * @return two Samples represented by the dataFile.
     */
    public static Sample[] getDataSet(String dataFile, String newick, AlphabetContainer con,
            int trainingSeqs, int testingSeqs, int splitseed)
            throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException,
            WrongSequenceTypeException {
		PhyloSample sample = SampleUtil.getDataSet(dataFile, newick, con);
		PhyloSample sampleTrain = sample.getSubSample(trainingSeqs, splitseed);
		PhyloSample sampleTest = PhyloSample.diff(sample, sampleTrain);
        if (testingSeqs > 0 && testingSeqs < sampleTest.getNumberOfElements()) {
			sampleTest = sampleTest.getSubSample(testingSeqs, splitseed);
        }
        return new Sample[] { sampleTrain, sampleTest };
    }

    /**
     * Loads four data sets from the given filenames. Converts those data sets to {@link PhyloSample}s. Returns the four
     * data sets in an array with: 0 => training foreground; 1 => testing foreground; 2 => training background;
     * 3 => testing background.
     * 
     * @param dataFileFGtrain
     * @param dataFileFGtest
     * @param dataFileBGtrain
     * @param dataFileBGtest
     * @param con
     * @param newick
     * @return array containing four samples.
     * @throws EmptySampleException
     * @throws WrongAlphabetException
     * @throws WrongLengthException
     * @throws IOException
     * @throws WrongSequenceTypeException
     */
    public static Sample[] getDataSets(
            String dataFileFGtrain, String dataFileFGtest, String dataFileBGtrain, String dataFileBGtest,
            AlphabetContainer con, String newick) throws EmptySampleException, WrongAlphabetException,
            WrongLengthException, IOException, WrongSequenceTypeException {
		LOGGER.info("Loading without splitting from fg.train, fg.test, bg.train, bg.test file");
        ArrayList<String> speciesOrder = Util.getOrderedArrayListFromTree(newick);
        SplitSequenceAnnotationParser parser = new SplitSequenceAnnotationParser("=", ";");

        Sample rawDataFGtrain = new Sample(con, new SparseStringExtractor(dataFileFGtrain, '>', parser));
		PhyloSample sampleTrainFG = new PhyloSample(rawDataFGtrain, speciesOrder);

        Sample rawDataFGtest = new Sample(con, new SparseStringExtractor(dataFileFGtest, '>', parser));
		PhyloSample sampleTestFG = new PhyloSample(rawDataFGtest, speciesOrder);

        Sample rawDataBGtrain = new Sample(con, new SparseStringExtractor(dataFileBGtrain, '>', parser));
		PhyloSample sampleTrainBG = new PhyloSample(rawDataBGtrain, speciesOrder);

        Sample rawDataBGtest = new Sample(con, new SparseStringExtractor(dataFileBGtest, '>', parser));
		PhyloSample sampleTestBG = new PhyloSample(rawDataBGtest, speciesOrder);

		return new PhyloSample[] { sampleTrainFG, sampleTestFG, sampleTrainBG, sampleTestBG };
    }

    /**
     * Loads two datasets from the given two filenames and converts them to {@link PhyloSample}s. Splits the two data
     * sets each into a training and a test dataset. Returns the four resulting datasets in an array with: 0 => training
     * foreground; 1 => testing foreground; 2 => training background; 3 => testing background.
     * 
     * @param dataFileFG
     * @param dataFileBG
     * @param con
     * @param newick
     * @param FILTER_TOP_SEQUENCES_BY_ANNOTATION
     * @param NUMBER_OF_SEQS_FOR_TRAIN_FG
     * @param NUMBER_OF_SEQS_FOR_TEST_FG
     * @param NUMBER_OF_SEQS_FOR_TRAIN_BG
     * @param NUMBER_OF_SEQS_FOR_TEST_BG
     * @param SPLIT_SEED
     * @return array containing four samples.
     * @throws FileNotFoundException
     * @throws WrongAlphabetException
     * @throws EmptySampleException
     * @throws WrongLengthException
     * @throws IOException
     * @throws WrongSequenceTypeException
     */
    @Deprecated
    public static Sample[] getDataSets(String dataFileFG, String dataFileBG, AlphabetContainer con, String newick,
            final int FILTER_TOP_SEQUENCES_BY_ANNOTATION,
            final int NUMBER_OF_SEQS_FOR_TRAIN_FG,
            final int NUMBER_OF_SEQS_FOR_TEST_FG,
            final int NUMBER_OF_SEQS_FOR_TRAIN_BG,
            final int NUMBER_OF_SEQS_FOR_TEST_BG,
            final int SPLIT_SEED)
            throws FileNotFoundException, WrongAlphabetException, EmptySampleException, WrongLengthException,
            IOException, WrongSequenceTypeException {
        System.out.println("Loading and splitting from a fg and bg file");
        ArrayList<String> speciesOrder = Util.getOrderedArrayListFromTree(newick);
        SplitSequenceAnnotationParser parser = new SplitSequenceAnnotationParser("=", ";");
        Sample rawDataFG = new Sample(con, new SparseStringExtractor(dataFileFG, '>', parser));
		PhyloSample sampleFG = new PhyloSample(rawDataFG, speciesOrder);
        if (FILTER_TOP_SEQUENCES_BY_ANNOTATION > 0) {
			sampleFG = sampleFG.getSubSampleByAnnotatedScore(0, FILTER_TOP_SEQUENCES_BY_ANNOTATION);
        }

		PhyloSample sampleTrainFG = sampleFG.getSubSample(NUMBER_OF_SEQS_FOR_TRAIN_FG, SPLIT_SEED);

		PhyloSample sampleTestFG = PhyloSample.diff(sampleFG, sampleTrainFG);
        if (NUMBER_OF_SEQS_FOR_TEST_FG > 0) {
			sampleTestFG = sampleTestFG.getSubSample(NUMBER_OF_SEQS_FOR_TEST_FG, SPLIT_SEED);
        }

        Sample rawDataBG = new Sample(con, new SparseStringExtractor(dataFileBG, '>', parser));
		PhyloSample sampleBG = new PhyloSample(rawDataBG, speciesOrder);
		PhyloSample sampleTrainBG = sampleBG.getSubSample(NUMBER_OF_SEQS_FOR_TRAIN_BG, SPLIT_SEED);
		PhyloSample sampleTestBG = PhyloSample.diff(sampleBG, sampleTrainBG);
        if (NUMBER_OF_SEQS_FOR_TEST_BG > 0) {
			sampleTestBG = sampleTestBG.getSubSample(NUMBER_OF_SEQS_FOR_TEST_BG, SPLIT_SEED);
        }
        return new Sample[] { sampleTrainFG, sampleTestFG, sampleTrainBG, sampleTestBG };
    }

    /**
     * Generates a fasta string from the given Sample.
     * 
     * @param sample
     * @return the corresponding fasta string
     */
    public static String toFasta(Sample sample) {
        StringBuilder sb = new StringBuilder();
        Sequence<?>[] seqs = sample.getAllElements();
        for (Sequence<?> seq : seqs) {
            if (seq instanceof MultiDimensionalDiscreteSequence) {
                for (int o = 0; o < ((MultiDimensionalDiscreteSequence) seq).getNumberOfSequences(); o++) {
                    Sequence<?> subseq = ((MultiDimensionalDiscreteSequence) seq).getSequence(o);
                    sb.append(">" + getAnnotationString(subseq.getAnnotation()) + "\n");
                    sb.append(subseq.toString() + "\n");
                }
            } else if (seq instanceof SimpleDiscreteSequence) {
                sb.append(">");
                sb.append(seq.toString());
            } else {
                throw new UnsupportedOperationException("Fasta of a PhyloSample not possible on sequences of type "
                        + seq.getClass());
            }
        }
        return sb.toString();
    }

    /**
     * Creates an annotation string from the given annotation array. The annotation delimiter is ";". The assign symbol
     * is "=".
     * 
     * @param annotation
     * @return
     */
    private static String getAnnotationString(SequenceAnnotation[] annotation) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < annotation.length; i++) {
            if (i > 0) {
                sb.append("; ");
            }
            sb.append(annotation[i].getType() + "=" + annotation[i].getIdentifier());
        }
        return sb.toString();
    }

    public static Sample[] getDataSets(Properties props, AlphabetContainer alphabet)
            throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException,
            WrongSequenceTypeException {
        String dataFileFG = Config.getProperty(props, "input.dataset.fg", false).asString();
        String dataFileBG = Config.getProperty(props, "input.dataset.bg", false).asString();
        String newick = Config.getPropertyWithFallBack(props, false, "model.bg.newick", "model.fg.newick",
                "model.newick").asString();
        int splitseed = Config.getProperty(props, "input.splitseed", new Random().nextInt() + "").asInt();

        Sample[] samples = new Sample[4];

        int train = Config.getProperty(props, "input.dataset.fg.train", "200").asInt();
        int test = Config.getProperty(props, "input.dataset.fg.test", "0").asInt();
        Sample[] tmp = SampleUtil.getDataSet(dataFileFG, newick, alphabet, train, test, splitseed);
        samples[0] = tmp[0];
        samples[1] = tmp[1];

        train = Config.getProperty(props, "input.dataset.bg.train", "500").asInt();
        test = Config.getProperty(props, "input.dataset.bg.test", "0").asInt();
        tmp = SampleUtil.getDataSet(dataFileBG, newick, alphabet, train, test, splitseed);
        
		samples[2] = tmp[0];
		samples[3] = tmp[1];

        return samples;
    }

	public static int[] getTruePositionSet(Sample sample, String posIdentifier) throws IOException {
		int[] truePos = new int[sample.getNumberOfElements()];
		Arrays.fill(truePos, -1);
		for (int i = 0; i < sample.getNumberOfElements(); i++) {
			SequenceAnnotation sequenceAnnotation = ((MultiDimensionalDiscreteSequence) sample.getElementAt(i)).getSequence(0).getAnnotation()[2];
			if (!sequenceAnnotation.getType().equals(posIdentifier)) {
				throw new IOException("Expected motif " + posIdentifier + ". But seen " + sequenceAnnotation.getType());
			}
			int pos = Integer.valueOf(sequenceAnnotation.getIdentifier());
			if (pos >= 0) {
				truePos[i] = pos;
			}
		}
		return truePos;
	}
}
