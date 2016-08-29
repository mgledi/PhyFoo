package io;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import de.jstacs.data.Sample;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.SubSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import util.AssociativeSort;

/**
 * This class represents a DataModel to extract an "important" subset from an
 * array of Sequences and a corresponding array of weights. It respects that
 * several sequences sharing one common parent. So each parent sequence can
 * contribute equal, independent of the number of subsequences from this parent.
 * 
 * @author Martin Nettling
 * 
 */
public class SequenceSpecificDataSelector {
    public Sequence<?>[][] seqsPerParent;
    public double weightProfile[][];
    public double normWeightProfile[][];

    public double[] weightsForTraining;
    @SuppressWarnings("rawtypes")
    public Sequence[] seqsForTraining;

    public SequenceSpecificDataSelector(Sample data, double[] weights) {
        transform(data, weights);
    }

    public void transform(Sample data, double[] weights) {
        // Sequence[] seqs = data.getAllElements();

        HashMap<Sequence<?>, ArrayList<Sequence<?>>> ss = new HashMap<Sequence<?>, ArrayList<Sequence<?>>>();
        HashMap<Sequence<?>, DoubleList> ws = new HashMap<Sequence<?>, DoubleList>();
        int i = 0;
        for (Sequence<?> s : data.getAllElements()) {
            Sequence<?> parent;
            if (s instanceof MultiDimensionalDiscreteSequence) {
                parent = s;
            } else {
                parent = ((SubSequence<?>) s).getParent();
            }
            if (!ss.containsKey(parent)) {
                ss.put(parent, new ArrayList<Sequence<?>>());
                ws.put(parent, new DoubleList());
            }
            ss.get(parent).add(s);

            ws.get(parent).add(weights == null ? 1. : weights[i]);
            i++;
        }
        seqsPerParent = new Sequence[ss.size()][];
        weightProfile = new double[ss.size()][];
        normWeightProfile = new double[ss.size()][];
        i = 0;
        for (Sequence<?> parent : ss.keySet()) {
            seqsPerParent[i] = new Sequence[ss.get(parent).size()];
            weightProfile[i] = new double[ss.get(parent).size()];
            normWeightProfile[i] = new double[ss.get(parent).size()];
            for (int k = 0; k < seqsPerParent[i].length; k++) {
                seqsPerParent[i][k] = ss.get(parent).get(k);
                weightProfile[i][k] = ws.get(parent).get(k);
            }
            AssociativeSort.quickSort(weightProfile[i], seqsPerParent[i]);
            i++;
        }

        for (i = 0; i < weightProfile.length; i++) {
            normWeightProfile[i] = Arrays.copyOf(weightProfile[i], weightProfile[i].length);
            Normalisation.sumNormalisation(normWeightProfile[i]);
        }
    }

    public void prepareForTraining(double t, double q) {
        ArrayList<Sequence<?>> seqs = new ArrayList<Sequence<?>>();
        DoubleList weights = new DoubleList();
        for (int i = 0; i < normWeightProfile.length; i++) {
            double sum = 0;
            for (int k = normWeightProfile[i].length - 1; k >= 0; k--) {
                sum += normWeightProfile[i][k];
                if (weightProfile[i][k] >= q) {
                    seqs.add(seqsPerParent[i][k]);
                    weights.add(weightProfile[i][k]);
                }
                if (sum > t) {
                    break;
                }
            }
        }
        weightsForTraining = weights.toArray();
        seqsForTraining = seqs.toArray(new Sequence[0]);
    }
}
