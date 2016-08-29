package io;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import de.jstacs.NotTrainedException;
import de.jstacs.WrongAlphabetException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.EmptySampleException;
import de.jstacs.data.Sample;
import de.jstacs.data.WrongLengthException;
import de.jstacs.data.alphabets.DNAAlphabet;
import de.jstacs.data.alphabets.DiscreteAlphabet;
import de.jstacs.data.alphabets.DoubleSymbolException;
import de.jstacs.data.sequences.ByteSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.data.sequences.SimpleDiscreteSequence;
import de.jstacs.data.sequences.WrongSequenceTypeException;
import de.jstacs.data.sequences.annotation.SequenceAnnotation;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner;
import de.jstacs.models.discrete.inhomogeneous.parameters.BayesianNetworkModelParameterSet;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.ToolBox;
import models.PhyloBackground;
import util.AssociativeSort;
import util.MatrixLinearisation;

/**
 * A Sample used in Phylogenetic approaches like Phylogenetic Footprinting to represent a set of alignments.
 * 
 * @author mnettling
 *
 */
public class PhyloSample extends Sample {

	private static final String SPECIES_IDENT = "species";
	private static final String GENE_IDENT = "gene";

	private final List<String> _species;
	
	public PhyloSample(String annotation, MultiDimensionalDiscreteSequence[] sequences, List<String> species)
	        throws EmptySampleException, WrongAlphabetException {
		super(annotation, sequences);
		_species = species;
	}

	/**
	 * Creates a new {@link PhyloSample} from a the given sample. Alignments are grouped by the given gene identifier.
	 * 
	 * @throws IOException
	 * @throws WrongSequenceTypeException
	 */
	public PhyloSample(Sample s, List<String> species)
	        throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException, WrongSequenceTypeException {
		super("multi-dimensional sample of " + s.getAnnotation(), getMultiDimensionalDiscreteSequences(s, SPECIES_IDENT, GENE_IDENT, species));
		_species = species;
	}

	/**
	 * generates a Set of Multidimensional Sequences from a "normal" Sample. A speciesIdentifier and a geneIdentifier must been set properly.
	 * 
	 * @throws IOException
	 * @throws WrongSequenceTypeException
	 */
	private static MultiDimensionalDiscreteSequence[] getMultiDimensionalDiscreteSequences(
			Sample data, 
			String speciesIdent, 
			String geneIdent,
	        List<String> speciesOrder) throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException, WrongSequenceTypeException {

		Hashtable<String, HashSet<String>> annot = data.getAnnotationTypesAndIdentifier();
		Hashtable<String, Integer> multiDimHashtable = ToolBox.parseHashSet2IndexHashtable(annot.get(speciesIdent));

		// check consistency if speciesOrder is set
		for (String species : speciesOrder) {
			if (!multiDimHashtable.containsKey(species)) {
				throw new IOException("A species-name (" + species + ") exisiting in the tree does not exist in the data-set.");
			}
		}

		// Reorders the hashmap, so that the order of the leafs in tree is exactly the same
		// as the order of sequences in the dataset
		int index = 0;

		for (String species : speciesOrder) {
			multiDimHashtable.put(species, index++);
		}

		// if there are other elements in the data set, which are not in the tree
		// set their indices for correctness
		for (Entry<String, Integer> entry : multiDimHashtable.entrySet()) {
			if (!speciesOrder.contains(entry.getKey())) {
				entry.setValue(index++);
			}
		}
		Hashtable<String, Integer> separateHashtable = ToolBox.parseHashSet2IndexHashtable(annot.get(geneIdent));
		int[][] matrix = data.getSequenceAnnotationIndexMatrix(geneIdent, separateHashtable, speciesIdent, multiDimHashtable);
		// SimpleDiscreteSequence[] content = new SimpleDiscreteSequence[multiDimHashtable.size()];
		SimpleDiscreteSequence[] content = new SimpleDiscreteSequence[speciesOrder.size()];
		MultiDimensionalDiscreteSequence[] seqs = new MultiDimensionalDiscreteSequence[matrix.length];
		SequenceAnnotation[] seqAn = new SequenceAnnotation[3];
		String id = null;
		String ref_id = null;
		for (int r = 0; r < matrix.length; r++) {
			for (int m = 0; m < speciesOrder.size(); m++) {
				// for (int m = 0; m < matrix[r].length; m++) {
				if (matrix[r][m] < 0) {
					throw new IllegalArgumentException();
				} else {
					content[m] = (SimpleDiscreteSequence) data.getElementAt(matrix[r][m]);
					if (m == 0) {
						id = "(" + content[m].getSequenceAnnotationByType(speciesIdent, 0).getIdentifier();
						ref_id = "(" + matrix[r][m];
					} else {
						id += ", " + content[m].getSequenceAnnotationByType(speciesIdent, 0).getIdentifier();
						ref_id += ", " + matrix[r][m];
					}
				}
			}
			getSequenceWithoutAlignmentGaps(content);

			seqAn[0] = new SequenceAnnotation(geneIdent, content[0].getSequenceAnnotationByType(geneIdent, 0).getIdentifier());
			seqAn[1] = new SequenceAnnotation(speciesIdent, id + ")");
			seqAn[2] = new SequenceAnnotation("SampleReference", ref_id + ")");
			seqs[r] = new MultiDimensionalDiscreteSequence(seqAn, content);
		}
		return seqs;
	}

	public List<String> getSpecies() {
		for (int i = 0; i < getNumberOfElements(); i++) {
		}
		return null;
	}

	/** @return the number of organism */
	public int getNumberOfOrganism() {
		return _species.size();
	}

	/**
	 * Transforms the given PhyloSample to a simple Sample consisting only of ByteSequences using the given Alphabet.
	 * 
	 * @param sample
	 * @param expectedAlphabet
	 * @param filter
	 * @return
	 * @throws WrongAlphabetException
	 * @throws WrongSequenceTypeException
	 * @throws EmptySampleException
	 */
	public static Sample toSingleSequenceDNASample(Sample sample, AlphabetContainer expectedAlphabet, SequenceAnnotation filter)
	        throws WrongAlphabetException, WrongSequenceTypeException, EmptySampleException {
		ArrayList<ByteSequence> seqs = new ArrayList<ByteSequence>();
		for (int i = 0; i < sample.getNumberOfElements(); i++) {
			if (!(sample.getElementAt(i) instanceof MultiDimensionalDiscreteSequence)) {
				throw new IllegalArgumentException("Can only convert " + MultiDimensionalDiscreteSequence.class);
			}
			MultiDimensionalDiscreteSequence seq = (MultiDimensionalDiscreteSequence) sample.getElementAt(i);
			// transform all MultiDimensionalDiscreteSequences to their gapless ByteSequeces
			for (int o = 0; o < seq.getNumberOfSequences(); o++) {
				SimpleDiscreteSequence[] tmp = new SimpleDiscreteSequence[] { (SimpleDiscreteSequence) seq.getSequence(o) };
				// remove gaps
				getSequenceWithoutAlignmentGaps(tmp);

				boolean filterMatch = false;
				if (filter != null) {
					for (SequenceAnnotation sa : tmp[0].getAnnotation()) {
						if (sa.getIdentifier().equals(filter.getIdentifier()) && sa.getType().equals(filter.getType())) {
							// this sequence matches the filter and should not be filtered
							filterMatch = true;
							break;
						}
					}
				}
				// generate correct ByteSequence (needed for BayesianNetworkModel
				tmp[0] = new ByteSequence(expectedAlphabet, tmp[0].getAnnotation(), tmp[0].toString(), "");
				if (filter == null || filterMatch) {
					seqs.add((ByteSequence) tmp[0]);
				}
			}
		}
		return new Sample("Transformed from PhyloSample" + sample.getAnnotation(), seqs.toArray(new ByteSequence[0]));
	}

	private static void getSequenceWithoutAlignmentGaps(SimpleDiscreteSequence[] content) throws WrongAlphabetException, WrongSequenceTypeException {
		String[] newSequences = new String[content.length];
		for (int o = 0; o < newSequences.length; o++) {
			newSequences[o] = new String();
		}

		DiscreteAlphabet alphabet = (DiscreteAlphabet) content[0].getAlphabetContainer().getAlphabetAt(0);
		for (int u = 0; u < content[0].getLength(); u++) {
			int sumgaps = 0;
			for (int o = 0; o < content.length; o++) {
				if (content[o].discreteVal(u) == alphabet.getCode("-")) {
					sumgaps++;
				}
			}

			// add if at least one symbol was seen
			if (sumgaps < content.length) {
				for (int o = 0; o < content.length; o++) {
					newSequences[o] += alphabet.getSymbolAt(content[o].discreteVal(u));
				}
			}
		}

		for (int o = 0; o < content.length; o++) {
			content[o] = new ByteSequence(content[o].getAlphabetContainer(), content[o].getAnnotation(), newSequences[o],
			        content[o].getAlphabetContainer().getDelim());
		}
	}

	/**
	 * This method parses the reference ids from the given String. Each sequence in a PhyloSample has at least one reference in an underlying sample
	 */
	public static int[] parseReferenceIds(String ref) {
		String[] refStrings = ref.substring(1, ref.length() - 1).split(",");
		int[] refIds = new int[refStrings.length];
		for (int i = 0; i < refIds.length; i++) {
			refIds[i] = Integer.valueOf(refStrings[i].trim());
		}
		return refIds;
	}

	/**
	 * The method trains for each {@link MultiDimensionalDiscreteSequence} in the given sample the branch lengths of given baseTopology using
	 * {@link PhyloBackground} of order 0. It stores the trained topology in the particular instance of {@link MultiDimensionalDiscreteSequence}.
	 * 
	 * @author Martin Nettling
	 */
	public static void train(PhyloSample sample, double weight, String baseTopology) throws Exception {
		PhyloBackground bg = new PhyloBackground((byte) 0, baseTopology, sample.getAlphabetContainer());
		PhyloBackground.ENABLE_CACHING = true;
		PhyloBackground.EDGE_LEARNING = true;

		double[] initialWeights = new double[sample.getNumberOfElements()];
		for (int i = 0; i < sample.getNumberOfElements(); i++) {
			initialWeights[i] = (1 - weight) / sample.getNumberOfElements();
		}

		for (int k = 0; k < 5; k++) {
			bg.train(sample, null);
		}
		String globalTopology = bg.getBNH().getVirtualTree(0).getNewickString();
		for (int i = 0; i < sample.getNumberOfElements(); i++) {
			System.out.println("Optimizing alignment " + (i + 1) + " of " + sample.getNumberOfElements());
			bg.reinitEdgelengths(baseTopology);
			double[] weights = Arrays.copyOf(initialWeights, initialWeights.length);
			weights[i] = weight;
			Sample trainingSample = new Sample("SubSample " + i, sample.getElementAt(i));
			bg.train(trainingSample, weights);
			String topology = combineBranchLengths(bg.getBNH().getVirtualTree(0).getNewickString(), globalTopology, weight);
			((MultiDimensionalDiscreteSequence) sample.getElementAt(i)).learnedTopology = topology;
		}

		PhyloBackground.ENABLE_CACHING = false;
		PhyloBackground.EDGE_LEARNING = false;
	}

	/**
	 * The method calculates the weighted mean for each branch length in the two given topologies.
	 * 
	 * @param newick1
	 * @param newick2
	 * @param weight
	 *            the weight for newick1
	 * @return
	 */
	public static String combineBranchLengths(String newick1, String newick2, double weight) {
		DoubleList dl1 = new DoubleList();
		Pattern p = Pattern.compile(":([0-9]*\\.[0-9]*)");
		Matcher m1 = p.matcher(newick1);
		while (m1.find()) {
			dl1.add(Double.valueOf(m1.group(1)));
		}

		DoubleList dl2 = new DoubleList();
		Matcher m2 = p.matcher(newick2);
		while (m2.find()) {
			dl2.add(Double.valueOf(m2.group(1)));
		}

		Matcher m3 = p.matcher(newick1);
		int last = 0;
		String combined = "";
		int s = 0;
		while (m3.find()) {
			for (int i = last; i < m3.start(1); i++) {
				combined += newick1.charAt(i);
			}
			combined += String.valueOf((dl1.get(s) * weight + dl2.get(s) * (1 - weight)));
			last = m3.end(1);
			s++;
		}
		for (int i = last; i < newick1.length(); i++) {
			combined += newick1.charAt(i);
		}
		return combined;
	}

	public static PhyloSample diff(PhyloSample data, PhyloSample... samples) throws EmptySampleException, WrongAlphabetException {
		Sample result = Sample.diff(data, samples);
		MultiDimensionalDiscreteSequence[] sequences = new MultiDimensionalDiscreteSequence[result.getNumberOfElements()];
		for (int i=0; i < result.getNumberOfElements(); i++) {
			if (result.getElementAt(i) instanceof MultiDimensionalDiscreteSequence) {
				sequences[i] = (MultiDimensionalDiscreteSequence) result.getElementAt(i);
			} else {
				throw new IllegalArgumentException("The resulting sample must only contain MultiDimensionalDiscreteSequences");
			}
		}
		return new PhyloSample(result.getAnnotation(), sequences, data._species);
	}

	/**
	 * 
	 * @param species_id
	 * @param numberOfSequencesToReturn
	 * @return
	 * @throws EmptySampleException
	 * @throws WrongAlphabetException
	 */
	public PhyloSample getSubSampleByAnnotatedScore(int species_id, int numberOfSequencesToReturn)
	        throws EmptySampleException, WrongAlphabetException {
		double[] maxScores = new double[getNumberOfElements()];
		// get all MultiDimensionalDiscreteSequences from PhyloSample
		MultiDimensionalDiscreteSequence[] seqs = new MultiDimensionalDiscreteSequence[getNumberOfElements()];
		for (int i = 0; i < seqs.length; i++) {
			seqs[i] = (MultiDimensionalDiscreteSequence) getElementAt(i);
			SequenceAnnotation[] annotation = seqs[i].getSequence(species_id).getAnnotation();
			boolean found = false;
			for (SequenceAnnotation a : annotation) {
				if (a.getType().equals("score")) {
					maxScores[i] = (Double.valueOf(a.getIdentifier()));
					found = true;
				}
			}
			if (!found) {
				throw new IllegalArgumentException("No annotation of type 'score' found.");
			}
		}

		// sort scores and sequences associated
		AssociativeSort.quickSort(maxScores, seqs);
		if (numberOfSequencesToReturn > seqs.length) {
			numberOfSequencesToReturn = seqs.length;
		}
		MultiDimensionalDiscreteSequence[] toReturn = Arrays.copyOfRange(seqs, seqs.length - numberOfSequencesToReturn, seqs.length);

		return new PhyloSample("Filtered Sample by Score Annotation ", toReturn, _species);
	}

	/**
	 * 
	 * @param species_id
	 * @param pwm
	 * @param sample
	 * @param percent
	 * @return
	 * @throws DoubleSymbolException
	 * @throws Exception
	 */
	public static Sample getSubSampleByPWMFilter(int species_id, double[][] pwm, PhyloSample sample, double percent, int numberOfSequencesToReturn)
	        throws DoubleSymbolException, Exception {
		BayesianNetworkModelParameterSet param = new BayesianNetworkModelParameterSet(new AlphabetContainer(new DNAAlphabet()), pwm.length, 0.0D,
		        "foreground model", StructureLearner.ModelType.IMM, (byte) 0, StructureLearner.LearningType.ML_OR_MAP);
		BayesianNetworkModel motif = new BayesianNetworkModel(param);

		motif.trained = true;
		int[][] structure = new int[pwm.length][1];
		for (int i = 0; i < structure.length; i++) {
			structure[i][0] = i;
		}
		motif.createConstraints(structure);
		for (int i = 0; i < pwm.length; i++) {
			motif.constraints[i].freq = pwm[i];
			motif.constraints[i].lnFreq = MatrixLinearisation.pi2lambda(pwm[i]);
		}

		double[] maxScores = new double[sample.getNumberOfElements()];
		// get all MultiDimensionalDiscreteSequences from PhyloSample
		MultiDimensionalDiscreteSequence[] seqs = new MultiDimensionalDiscreteSequence[sample.getNumberOfElements()];
		for (int i = 0; i < seqs.length; i++) {
			seqs[i] = (MultiDimensionalDiscreteSequence) sample.getElementAt(i);
		}

		for (int i = 0; i < seqs.length; i++) {
			// generate Array of size 1 (API of getSequenceWithoutAlignment)
			SimpleDiscreteSequence[] tmp = new SimpleDiscreteSequence[] { (SimpleDiscreteSequence) seqs[i].getSequence(species_id) };
			// remove gaps
			getSequenceWithoutAlignmentGaps(tmp);
			// generate correct ByteSequence (needed for BayesianNetworkModel
			tmp[0] = new ByteSequence(new AlphabetContainer(new DNAAlphabet()), tmp[0].toString());
			// get best PWM Score
			maxScores[i] = getMaximumPWMScore(tmp[0], motif);
		}

		// sort scores and sequences associated
		AssociativeSort.quickSort(maxScores, seqs);
		int size = (int) (sample.getNumberOfElements() * Math.abs(percent));
		if (numberOfSequencesToReturn == 0) {
			numberOfSequencesToReturn = size;
		}

		// generate set with all putative sequence IDs
		// if percent is smaller than 0, return "bad" sequences, else good sequences
		LinkedList<Integer> idSet = new LinkedList<Integer>();
		for (int i = 0; i < size; i++) {
			idSet.add(percent < 0 ? i : (seqs.length - 1 - i));
		}

		ArrayList<MultiDimensionalDiscreteSequence> randomSubSet = new ArrayList<MultiDimensionalDiscreteSequence>();
		Random r = new Random();
		for (int i = 0; i < numberOfSequencesToReturn && idSet.size() > 0; i++) {
			int element = r.nextInt(idSet.size());
			randomSubSet.add(seqs[idSet.get(element)]);
			idSet.remove(element);
		}
		return new Sample("Filtered Sample", randomSubSet.toArray(new MultiDimensionalDiscreteSequence[0]));
	}

	private static double getMaximumPWMScore(SimpleDiscreteSequence seq, BayesianNetworkModel bnm) throws NotTrainedException, Exception {
		double max = Double.NEGATIVE_INFINITY;
		for (int u = 0; u <= seq.getLength() - bnm.getLength(); u++) {
			double val = bnm.getLogProbFor(seq.getSubSequence(u, bnm.getLength()));
			if (val > max) {
				max = val;
			}
		}
		return max;
	}

	/**
	 * @param sequences
	 * @param s
	 * @param splitSeed
	 * 
	 * @throws WrongSequenceTypeException
	 * @throws IOException
	 * @throws WrongLengthException
	 * @throws WrongAlphabetException
	 * @throws EmptySampleException
	 * 
	 */
	public PhyloSample getSubSample(int sequences, int splitSeed)
	        throws EmptySampleException, WrongAlphabetException, WrongLengthException, IOException, WrongSequenceTypeException {
		// generate set with all sequence IDs
		LinkedList<Integer> idSet = new LinkedList<Integer>();
		for (int i = 0; i < getNumberOfElements(); i++) {
			idSet.add(i);
		}
		MultiDimensionalDiscreteSequence[] result = new MultiDimensionalDiscreteSequence[sequences];
		Random r = new Random(splitSeed);
		for (int i = 0; i < sequences && idSet.size() > 0; i++) {
			int element = r.nextInt(idSet.size());
			result[i] = (MultiDimensionalDiscreteSequence) getElementAt(idSet.get(element));
			idSet.remove(element);
		}
		return new PhyloSample("Subsample of " + getAnnotation(), result, _species);
	}

	public PhyloSample shuffle(String newick)
	        throws WrongAlphabetException, WrongSequenceTypeException, WrongLengthException, EmptySampleException, IOException {
		Random rand = new Random();

		MultiDimensionalDiscreteSequence[] result = new MultiDimensionalDiscreteSequence[getNumberOfElements()];
		for (int i = 0; i < getNumberOfElements(); i++) {
			MultiDimensionalDiscreteSequence source = (MultiDimensionalDiscreteSequence) getElementAt(i);
			DiscreteAlphabet alphabet = (DiscreteAlphabet) getAlphabetContainer().getAlphabetAt(0);
			int[] shuffle = new int[source.getLength()];
			for (int u = 0; u < shuffle.length; u++) {
				shuffle[u] = u;
			}
			// determine new order
			for (int u = 0; u < source.getLength(); u++) {
				int swap = rand.nextInt(source.getLength());
				int tmp = shuffle[u];
				shuffle[u] = shuffle[swap];
				shuffle[swap] = tmp;
			}

			SimpleDiscreteSequence[] singleSeq = new SimpleDiscreteSequence[_species.size()];
			for (int o = 0; o < source.getNumberOfSequences(); o++) {
				String tmpSeq = new String();
				for (int u : shuffle) {
					tmpSeq += alphabet.getSymbolAt(source.getSequence(o).discreteVal(u));
				}
				singleSeq[o] = new ByteSequence(getAlphabetContainer(), source.getSequence(o).getAnnotation(), tmpSeq, getAlphabetContainer().getDelim());
			}
			result[i] = new MultiDimensionalDiscreteSequence(source.getAnnotation(), singleSeq);
		}
		return new PhyloSample("Shuffled Sample of " + getAnnotation(), result,  _species);
	}
}
