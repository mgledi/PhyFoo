package models;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.List;
import java.util.Properties;

import org.apache.log4j.Logger;

import com.google.common.base.Splitter;

import de.jstacs.NonParsableException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.algorithms.optimization.termination.AbstractTerminationCondition;
import de.jstacs.algorithms.optimization.termination.SmallDifferenceOfFunctionEvaluationsCondition;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.io.FileManager;
import de.jstacs.io.XMLParser;
import de.jstacs.models.AbstractModel;
import de.jstacs.models.Model;
import de.jstacs.models.mixture.AbstractMixtureModel.Parameterization;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.models.mixture.motif.positionprior.UniformPositionPrior;
import de.jstacs.parameters.SimpleParameter.IllegalValueException;
import util.Config;

/**
 * This class provides methods to create and handle different types of models given some {@link Properties}. It helps to
 * configure and to transform them.
 * 
 * @author Martin Nettling
 * 
 */
public class ModelUtil {
	private static Logger LOGGER = Logger.getLogger(ModelUtil.class);

    /**
     * This method instantiates a new untrained {@link SingleHiddenMotifMixture} from the given properties. The
     * following properties are taken into account:
     * <ul>
     * <li>motif length
     * <li>order of motif and flanking model
     * <li>newick string for motif and flanking model
     * <li>initial zoops parameter
     * <li>usage of {@link StrandModel} for motif or not
     * <li>usage of {@link AbstractPhyloModel}s or {@link AbstractAlignmentBasedModel}
     * </ul>
     * 
     * @param props
     * @param alphabet
     * @param stopCondition
     * @return a new {@link SingleHiddenMotifMixture}
     * @throws IllegalArgumentException
     * @throws IllegalValueException
     * @throws CloneNotSupportedException
     * @throws WrongAlphabetException
     * @throws IOException
     * @throws NonParsableException
     */
    public static SingleHiddenMotifMixture getNewSHM(Properties props, AlphabetContainer alphabet,
            AbstractTerminationCondition stopCondition)
            throws IllegalArgumentException, IllegalValueException, CloneNotSupportedException, WrongAlphabetException,
            NonParsableException, IOException {
        if (Config.getProperty(props, "algorithm.continue", "false").asBoolean() &&
                new File(Config.getProperty(props, "input.model.shm", "").asString()).exists()) {
			LOGGER.info("Try to continue training on existing SHM");
            return getTrainedSHM(props);
        }

        int motifLength = Config.getProperty(props, "model.fg.length", false).asInt();
        byte motifOrder = Config.getProperty(props, "model.fg.order", "0").asByte();
        String motifNewick = Config.getPropertyWithFallBack(props, false, "model.fg.newick", "model.newick").asString();
        byte flankingOrder =
                Config.getPropertyWithFallBack(props, "0", "model.bg.flanking.order", "model.bg.order").asByte();
        String flankingNewick = Config.getPropertyWithFallBack(props, false, "model.bg.flanking.newick",
                "model.bg.newick", "model.newick").asString();
        boolean useStrandModel = Config.getProperty(props, "model.strandmodel", "false").asBoolean();
        boolean usePhyloModels = Config.getProperty(props, "model.evolution", "true").asBoolean();
        double zoops = Config.getProperty(props, "algorithm.zoops", "-1").asDouble();
        Model motif;
        Model flanking;
        if (!usePhyloModels) {
            motif = new AlignmentBasedModel(motifLength, motifOrder, alphabet);
            flanking = new AlignmentBasedBGModel(flankingOrder, alphabet);
        } else {
            motif = new PhyloBayesModel(motifLength, motifOrder, motifNewick, alphabet);
            flanking = new PhyloBackground(flankingOrder, flankingNewick, alphabet);
        }

		if (Config.getProperty(props, "model.fg.condprobs", true) != null) {
			double[][] condProbs = ((PhyloPreparedAbstractModel)motif).getCondProbs();
			parseCondProbs(Config.getProperty(props, "model.fg.condprobs", true).asString(), condProbs);
			((PhyloPreparedAbstractModel) motif).setCondProbs(condProbs);
		}

        if (useStrandModel) {
            motif = new StrandModel(
                    motif,
                    1, // we need only one iteration, because we use EMHandler
                    new double[] { 1, 1 },
                    1d,
                    new SmallDifferenceOfFunctionEvaluationsCondition(1E-4),
                    Parameterization.LAMBDA
                    );
            // the flanking model must not be encapsulated if all sequences are given from the same strand
        }
        SingleHiddenMotifMixture shm;
        if (zoops > 0 && zoops <= 1) {
            shm = new SingleHiddenMotifMixture(
                    motif,
                    flanking,
                    true,
                    1, // we need only one iteration, because we use EMHandler
                    zoops,
                    new UniformPositionPrior(),
                    1,
                    stopCondition,
                    Parameterization.LAMBDA
                    );
        } else {
            shm = new SingleHiddenMotifMixture(
                    motif,
                    flanking,
                    true,
                    1, // we need only one iteration, because we use EMHandler
                    new double[] { 1d, 1d },
                    new UniformPositionPrior(),
                    1,
                    stopCondition,
                    Parameterization.LAMBDA
                    );
        }


        shm.setShiftCorrection(Config.getProperty(props, "algorithm.enableshiftcorrection", "false").asBoolean());
        
        return shm;
    }

	private static void parseCondProbs(String probsString, double[][] condProbs) {
		List<String> elements = Splitter.on(",").trimResults().splitToList(probsString);
		Iterator<String> elementsIterator = elements.iterator();
		for (int i = 0; i < condProbs.length; i++) {
			for (int a = 0; a < condProbs[i].length; a++) {
				condProbs[i][a] = Double.valueOf(elementsIterator.next());
			}
		}
	}

	/**
	 * This method instantiates a new untrained background model from the given properties. The following properties are taken into account:
	 * <ul>
	 * <li>order
	 * <li>newick string
	 * <li>usage of {@link PhyloBackground}s or {@link AlignmentBasedBGModel}
	 * </ul>
	 * 
	 * @param props
	 * @param alphabet
	 * @return a new background model
	 * @throws IllegalArgumentException
	 * @throws IllegalValueException
	 * @throws CloneNotSupportedException
	 * @throws WrongAlphabetException
	 * @throws IOException
	 * @throws NonParsableException
	 */
    public static AbstractModel getNewBackground(Properties props, AlphabetContainer alphabet)
            throws IllegalArgumentException, IllegalValueException, CloneNotSupportedException, WrongAlphabetException,
            NonParsableException, IOException {
        if (Config.getProperty(props, "algorithm.continue", "false").asBoolean() &&
                new File(Config.getProperty(props, "input.model.bg", "").asString()).exists()) {
			LOGGER.info("Try to continue training on existing background");
            return getTrainedBackground(props);
        }
        byte order = Config.getProperty(props, "model.bg.order", "0").asByte();
        String newick = Config.getPropertyWithFallBack(props, false, "model.bg.newick", "model.newick").asString();
        boolean useStrandModel = Config.getProperty(props, "model.strandmodel", "false").asBoolean();
        boolean usePhyloModels = Config.getProperty(props, "model.evolution", "true").asBoolean();

        AbstractModel bg;
        if (!usePhyloModels) {
            bg = new AlignmentBasedBGModel(order, alphabet);
        } else {
            bg = new PhyloBackground(order, newick, alphabet);
        }
        if (useStrandModel) {
            // the flanking model must not be encapsulated if all sequences are given form the same strand
        }
        return bg;
    }

    /**
     * Accesses the parameter input.model.background and loads the corresponding {@link PhyloBackground} or
     * {@link AlignmentBasedBGModel}.
     * 
     * @return the already trained {@link PhyloPreparedAbstractModel} that is a {@link PhyloBackground} or
     *         {@link AlignmentBasedBGModel}
     * @throws NonParsableException
     * @throws IOException
     */
    public static PhyloPreparedAbstractModel getTrainedBackground(Properties props) throws NonParsableException,
            IOException {
        String trainedBGfile = Config.getProperty(props, "input.model.background", null).asString();
        StringBuffer xml = FileManager.readFile(new File(trainedBGfile));
        PhyloPreparedAbstractModel bgModel = null;
        // if a BayesNetHandler can be found this model is a PhlyoModel
        // this method is not able to handle a BG-model in a strandmodel
        if (XMLParser.hasTag(xml, "_bnh", null, null)) { // ugly hack
            bgModel = new PhyloBackground(xml);
        } else {
            bgModel = new AlignmentBasedBGModel(xml);
        }
        if (!bgModel.isTrained()) {
            throw new IllegalArgumentException("Tried to load an untrained BackgroundModel");
        }
        return bgModel;
    }

    /**
     * Accesses the parameter input.model.shm and loads the corresponding {@link SingleHiddenMotifMixture}.
     * 
     * @return the already trained {@link SingleHiddenMotifMixture}
     * @throws NonParsableException
     * @throws IOException
     */
    public static SingleHiddenMotifMixture getTrainedSHM(Properties props) throws NonParsableException, IOException {
        String trainedSHMfile = Config.getProperty(props, "input.model.shm", false).asString();
        StringBuffer xml = FileManager.readFile(new File(trainedSHMfile));
        SingleHiddenMotifMixture bestShm = new SingleHiddenMotifMixture(xml);
        if (!bestShm.isTrained()) {
            throw new IllegalArgumentException("Tried to load an untrained SingleHiddenMotifMixture");
        }
        return bestShm;
    }
}
