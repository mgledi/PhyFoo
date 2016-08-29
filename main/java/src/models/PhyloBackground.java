package models;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.apache.log4j.Logger;

import algorithm.MeanFieldForBayesNet;
import algorithm.SimpleNodeElimination;
import bayesNet.BayesNet;
import bayesNet.BayesNetHandler;
import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.io.XMLParser;
import de.jstacs.results.NumericalResultSet;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import evolution.EvolModel;
import io.Alphabet;
import optimizing.AbstractFreeEnergyOptimizer;
import optimizing.FE_byParameter_ForBayesNodeJstacs;
import optimizing.PreparedConditions;
import optimizing.TemperatureOptimizer;
import optimizing.branchlengths.EdgeOptimizer;
import util.Util;

/**
 *
 **/
public class PhyloBackground extends AbstractPhyloModel {
	private static Logger LOGGER = Logger.getLogger(PhyloBackground.class);

	/** If filtered data is assumed another optimizer must be used */
	public static int ASSUME_FILTERED_DATA = 0;

	public static boolean TRAIN_MOTIF = true;
	/** Shows if the edgle-lenghts should be learned during training */
	public static boolean EDGE_LEARNING = false;
	public static boolean CALC_LOGLIKELIHOOD = false;
	public static boolean CLONING_ALLOWED = true;

	/** */
	private boolean trained = true;

	/**
	 * Constructor
	 * 
	 * @param order
	 * @param newickString
	 * 
	 */
	public PhyloBackground(byte order, String newickString, AlphabetContainer con) {
		super(con, order, 0);
		bnh = buildStructure(newickString);
		fillParameter(false);
	}

	public PhyloBackground(StringBuffer xml) throws NonParsableException {
		super(xml);
	}

	public void preTrain(Sample data) {
		int org = ((MultiDimensionalDiscreteSequence) data.getElementAt(0)).getEmptyContainer().length;
		int a = (int) data.getAlphabetContainer().getAlphabetAt(0).length() - 1;
		a = 4; // TODO: variable
		for (int ord = 0; ord <= order; ord++) {
			int[][] con = new int[ord + 1][org];
			int[] subcon = new int[ord + 1];
			double[] counts = new double[(int) Math.pow(a, ord + 1)];
			int hash;
			for (int i = 0; i < data.getNumberOfElements(); i++) {
				for (int u = 0; u < data.getElementAt(i).getLength() - ord; u++) {
					for (int k = 0; k < ord + 1; k++) {
						((MultiDimensionalDiscreteSequence) data.getElementAt(i)).fillContainer(con[k], u + k);
					}

					for (int o = 0; o < org; o++) {
						hash = 0;
						for (int k = 0; k < ord + 1; k++) {
							subcon[k] = con[k][o];
							if (subcon[k] == 4) {
								hash = -1;
								break;
							}
							hash += subcon[k] * Math.pow(a, ord - k);
						}
						if (hash > -1) {
							counts[hash]++;
						}
					}
				}
			}
			for (int i = 0; i < counts.length; i += a) {
				double[] sub = Arrays.copyOfRange(counts, i, i + a);
				Normalisation.sumNormalisation(sub);
				for (int k = 0; k < a; k++) {
					counts[i + k] = sub[k];
				}
			}
			bnh.getVirtualTree(ord).getEvolModel().setStatDistr(counts);
			bnh.getVirtualTree(ord).initParametersFromGF();
		}
	}

	@Override
	public void train(Sample data, double[] weights) throws Exception {
		trained = true;
		long startTime = System.currentTimeMillis();
		// extends the weights from sequencelevel to nucleotide level
		double[] internalWeigths = weights != null ? extendWeights(data, weights) : null;
		int size = 0;
		for (int i = 0; i < data.getNumberOfElements(); i++) {
			size += data.getElementAt(i).getLength() - order;
		}

		/* singlethreading */
		MeanFieldForBayesNet[] tmpPointer = new MeanFieldForBayesNet[size];
		// zerlege jede Sequenz des samples in alle Teilsequenzen
		Sequence<?> subseq;
		int p = 0;
		LOGGER.info("Optimizing " + tmpPointer.length + " MeanFields of background.");
		for (int i = 0; i < data.getNumberOfElements(); i++) {
			for (int u = 0; u < data.getElementAt(i).getLength() - order; u++) {
				subseq = data.getElementAt(i).getSubSequence(u, order + 1);
				tmpPointer[p] = getMFfromCache(subseq);
				tmpPointer[p].optimizeByNormalisation();
				p++;
			}
		}
		LOGGER.info("Optimized " + tmpPointer.length + " MeanFieldForBayesNet for background.");

		if (TRAIN_MOTIF) {
			LOGGER.info("Use AbstractFreeEnergyOptimizer.");
			AbstractFreeEnergyOptimizer optPar = new FE_byParameter_ForBayesNodeJstacs();
			optPar.setModel(this);
			optPar.setVariationalLikelihood(tmpPointer);
			optPar.setTerminationCondition(PreparedConditions.SMALL_DIFFERENCE_OF_FUNCTIONS.clone());
			optPar.setWeights(internalWeigths);

			for (int k = 0; k < bnh.motifLength; k++) {
				((FE_byParameter_ForBayesNodeJstacs) optPar).initPosition(k);
				optPar.startOptimizing(false);
			}

			if (ASSUME_FILTERED_DATA > 0) {
				LOGGER.info("Use TemperatureOptimizer.");
				optPar = new TemperatureOptimizer(false);
				optPar.setVariationalLikelihood(tmpPointer);
				optPar.setModel(this);
				optPar.startOptimizing(false);
				LOGGER.info("Learned Temperture (0) = " + bnh.getVirtualTree(0).TEMPERATURE);
			}
		}

		if (EDGE_LEARNING) {
			EdgeOptimizer edgeOptimizer = new EdgeOptimizer(false, EDGE_LEARNING_ONE_LENGTH);
			edgeOptimizer.setModel(this);
			edgeOptimizer.setVariationalLikelihood(tmpPointer);
			edgeOptimizer.setWeights(internalWeigths);
			for (int k = 0; k < bnh.motifLength; k++) {
				edgeOptimizer.initPosition(k);
				edgeOptimizer.startOptimizing(false);
			}

			LOGGER.info("----- New NewickStrings PhyloBackground (" + (order + 1) + ")");
			for (int k = 0; k <= order; k++) {
				LOGGER.info(bnh.getVirtualTree(k).getNewickString());
			}
		}
		LOGGER.info("Finished Training in " + (System.currentTimeMillis() - startTime) + " milli seconds.");
	}

	private double[] extendWeights(Sample data, double[] weights) {
		DoubleList extended = new DoubleList();
		for (int i = 0; i < data.getNumberOfElements(); i++) {
			for (int k = 0; k < data.getElementAt(i).getLength(); k++) {
				extended.add(weights[i]);
			}
		}
		return extended.toArray();
	}

	public double getLogProbFor(Sequence sequence, int startpos, int endpos) throws IllegalArgumentException, NotTrainedException {
		return getLogProbFor((MultiDimensionalDiscreteSequence) sequence, startpos, endpos);
	}

	public byte getMaximalMarkovOrder() {
		return (byte) (order);
	}

	// this caching speedups the Calculation of getLogProbFor by factor of five
	private HashMap<MultiDimensionalDiscreteSequence, HashMap<Integer, HashMap<Integer, Double>>> _cachedScores = new HashMap<MultiDimensionalDiscreteSequence, HashMap<Integer, HashMap<Integer, Double>>>();

	/** cleans the cache */
	public void cleanCache() {
		_cachedScores.clear();
	}

	/** true: if a score was already computed */
	public boolean isScoreCached(MultiDimensionalDiscreteSequence sequence, int startpos, int endpos) {
		if (_cachedScores.containsKey(sequence) && _cachedScores.get(sequence).containsKey(startpos)
		        && _cachedScores.get(sequence).get(startpos).containsKey(endpos)) {
			return true;
		} else {
			return false;
		}
	}

	/** inits needed HashMaps for storing a score */
	public void prepareCache(MultiDimensionalDiscreteSequence sequence, int startpos, int endpos, double score) {
		if (!_cachedScores.containsKey(sequence)) {
			_cachedScores.put(sequence, new HashMap<Integer, HashMap<Integer, Double>>());
		}

		if (!_cachedScores.get(sequence).containsKey(startpos)) {
			_cachedScores.get(sequence).put(startpos, new HashMap<Integer, Double>());
		}

		_cachedScores.get(sequence).get(startpos).put(endpos, score);
	}

	public long count = 0;
	public static long calls = 0;

	/**
	 * Calculates the log probability for the given {@link MultiDimensionalDiscreteSequence} from starpos to endpos
	 * 
	 * @param sequence
	 * @param startpos
	 * @param endpos
	 * @return
	 * @throws NotTrainedException
	 * @throws IllegalArgumentException
	 * @throws Exception
	 */
	public double getLogProbFor(MultiDimensionalDiscreteSequence sequence, int startpos, int endpos) throws IllegalArgumentException, NotTrainedException {
		if (endpos < startpos) {
			return 0;
		}
		this.check(sequence, startpos, endpos);
		// if the score for this sequence and this start and endpos was already calculated, return it.
		if (isScoreCached(sequence, startpos, endpos)) {
			return _cachedScores.get(sequence).get(startpos).get(endpos);
		}
		String oldToplogy = this.getBNH().getVirtualTree(0).getNewickString();
		if (sequence.learnedTopology != null) {
			this.reinitEdgelengths(sequence.learnedTopology);
		}

		Sequence<?> subseq = null;
		int seqLength = endpos - startpos + 1; // startpos and endpos are inclusive

		ArrayList<MeanFieldForBayesNet> tmpPointer = new ArrayList<MeanFieldForBayesNet>();
		// first position
		subseq = sequence.getSubSequence(startpos, Math.min(order + 1, seqLength));
		tmpPointer.add(getMFfromCache(subseq));
		// divide the sequence in subsequences of length of order
		for (int u = startpos + 1; u <= endpos - order; u++) {
			subseq = sequence.getSubSequence(u, order + 1);
			tmpPointer.add(getMFfromCache(subseq));
		}

		// calculate log probability
		double logSum = 0;
		tmpPointer.get(0).initObservation();
		tmpPointer.get(0).init();
		tmpPointer.get(0).optimizeByNormalisation();
		for (int i = 0; i < order && startpos <= endpos; i++, startpos++) {
			if (CALC_LOGLIKELIHOOD) {
				SimpleNodeElimination.instantiate(bnh.getNet());
				logSum -= tmpPointer.get(0).calcLogLikelihood();
				// logSum -= tmpPointer.get(0).calcLogLikelihood(i-1, i); //TODO: NodeSumming is buggy
			} else {
				logSum += tmpPointer.get(0).calcFreeEnergy(i, i);
			}
			count++;
			// System.out.println("1 " + " " + pref + " = " + -logSum);
		}
		for (int i = 0; i < tmpPointer.size() && startpos <= endpos; i++, startpos++) {
			tmpPointer.get(i).initObservation();
			tmpPointer.get(i).init();
			tmpPointer.get(i).optimizeByNormalisation();
			if (CALC_LOGLIKELIHOOD) {
				SimpleNodeElimination.instantiate(bnh.getNet());
				logSum -= tmpPointer.get(i).calcLogLikelihood();
				// logSum -= tmpPointer.get(i).calcLogLikelihood(); //TODO: NodeSumming is buggy
			} else {
				logSum += tmpPointer.get(i).calcFreeEnergy(order, order);
			}
			count++;
			// System.out.println("2 " + " " + pref + " = " + -logSum);
		}
		prepareCache(sequence, startpos, endpos, -logSum);
		// System.out.println(pref + " = " + -logSum);

		if (sequence.learnedTopology != null) {
			this.reinitEdgelengths(oldToplogy);
		}
		return -logSum;
	}

	@Override
	public double getLogPriorTerm() throws Exception {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public String getInstanceName() {
		return "PhyloBackground (" + order + ")";
	}

	/** true, if the model is trained (filled with values */
	public boolean isTrained() {
		return trained;
	}

	@Override
	public NumericalResultSet getNumericalCharacteristics() throws Exception {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public StringBuffer toXML() {
		StringBuffer xml = bnh.toXML();
		XMLParser.addTags(xml, "_bnh");
		XMLParser.appendObjectWithTags(xml, order, "order");
		XMLParser.appendObjectWithTags(xml, alphabets, "alphabets");
		return xml;
	}

	@Override
	protected void fromXML(StringBuffer xml) throws NonParsableException {
		length = 0;
		alphabets = (AlphabetContainer) XMLParser.extractObjectForTags(xml, "alphabets");
		order = (Byte) XMLParser.extractObjectForTags(xml, "order");
		bnh = new BayesNetHandler(XMLParser.extractForTag(xml, "_bnh"));
		bnh.setSimpleRoot();
		trained = true;
	}

	@Override
	public PhyloBackground clone() {
		if (!CLONING_ALLOWED)
			return this;
		StringBuffer xml = this.toXML();
		try {
			return new PhyloBackground(xml);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return this;
	}

	// #########################################################################
	// ########################### Hilfsfunktionen #############################
	// #########################################################################
	public void fillParameter(boolean random) {
		EvolModel gfs = EvolModel.getInstance(1);
		for (int i = 0; i < bnh.getVirtualTrees().length; i++) {
			gfs = EvolModel.getInstance((int) Math.pow(Alphabet.size, bnh.getVirtualTree(i).getNode(0).numberOfParents));
			bnh.getVirtualTree(i).setEvolModel(gfs);
			if (random) {
				bnh.getVirtualTree(i).getEvolModel().setStatDistr(Util.getRandomStochMatrix(gfs.getDimension(), Alphabet.size));
			} else {
				bnh.getVirtualTree(i).getEvolModel().setStatDistr(Util.getEqualStochMatrix(gfs.getDimension(), Alphabet.size));
			}
			bnh.getVirtualTree(i).initParametersFromGF();
		}
	}

	/**
	 * baut die Graphische Struktur f�r ein Motiv mit �bergebener l�nge zu dem phyloBaum repr�sentiert durch newick
	 */
	private BayesNetHandler buildStructure(String newick) {
		// Struktur des Neztes bauen <- ACHTUNG: NOCH KEINE DATEN
		BayesNet wholeNet = new BayesNet("BG");
		BayesNetHandler bnh = new BayesNetHandler(wholeNet);

		for (int i = 0; i <= order; i++) {
			bnh.addBayesNet(io.NewickToBayesNet.getTree(newick, "bg_" + i), newick);
			// mit Vorg�nger verkn�pfen
			for (int j = 0; j < order; j++) {
				if (i - j > 0) {
					bnh.connectVirtualTrees(bnh.getVirtualTree(i - j - 1), bnh.getVirtualTree(i));
				}
			}
		}

		bnh.setSimpleRoot();
		return bnh;
	}

	/** Emitiert ein zuf�lliges Alignment */
	public String[] emitSample(int length) {
		// TODO: iMM(1+)
		int org = bnh.getVirtualTree(0).numberOfLeafs;
		String[] seq = new String[org];
		for (int o = 0; o < org; o++) {
			seq[o] = "";
		}
		// generiere eine initiale vollst�ndige Beobachtung f�r das Netz
		int[] oldObs;
		int[] fullObs = bnh.drawFullObservation();
		for (int i = 0; i <= order && i < length; i++) {
			for (int o = 0; o < org; o++) {
				alphabets.getSymbol(0, fullObs[bnh.getVirtualTree(i).getLeaf(o).nodeNumber]);
			}
		}

		// generiere f�r alle weiteren Positionen die Volls�ndigen beobachtungen
		// mit Erinnerung an die Alte
		for (int u = order + 1; u < length; u++) {
			oldObs = new int[bnh.getNet().numberOfNodes];
			Arrays.fill(oldObs, -1);
			// verschieben der Beobachtung um eine Position nach links
			for (int i = 1; i <= order; i++) {
				for (int k = 0; k < bnh.getVirtualTree(i).numberOfNodes; k++) {
					oldObs[bnh.getVirtualTree(i - 1).getNode(k).nodeNumber] = fullObs[bnh.getVirtualTree(i).getNode(k).nodeNumber];
				}
			}
			for (int o = 0; o < org; o++) {
				alphabets.getSymbol(0, fullObs[bnh.getVirtualTree(order).getLeaf(o).nodeNumber]);
			}
			fullObs = bnh.drawFullObservation(oldObs);
		}
		return seq;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.jstacs.models.discrete.DiscreteGraphicalModel#check(de.jstacs.data .Sequence, int, int)
	 */
	protected void check(Sequence<?> sequence, int startpos, int endpos) throws NotTrainedException, IllegalArgumentException {
		if (!trained) {
			throw new NotTrainedException();
		} else if (!alphabets.checkConsistency(sequence.getAlphabetContainer().getSubContainer(startpos, endpos - startpos + 1))) {
			throw new IllegalArgumentException("This sequence is not possible with the given alphabet.");
		} else if (startpos < 0) {
			throw new IllegalArgumentException("This startposition is impossible. Try: 0 <= startposition");
		} else if (startpos > endpos || endpos >= sequence.getLength()) {
			throw new IllegalArgumentException("This endposition is impossible. Try: startposition(" + startpos + ") <= endposition(" + endpos
			        + ") < sequence.length (" + sequence.getLength() + ")");
		} else if (endpos - startpos + 1 != length && length > 0) {
			throw new IllegalArgumentException("This sequence has not length " + length + ".");
		}
	}
}
