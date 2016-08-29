package classification;

import java.util.Arrays;
import java.util.HashSet;

import de.jstacs.classifier.MeasureParameters;
import de.jstacs.classifier.modelBased.ModelBasedClassifier;
import de.jstacs.data.Sample;
import de.jstacs.models.AbstractModel;
import de.jstacs.results.Result;
import de.jstacs.results.ResultSet;
import de.jstacs.utils.IntList;
import io.PhyloSample;
import util.IntPair;

public class ClassificationUtil {

    public static Result[] performClassification(AbstractModel shm, AbstractModel myBG,
            Sample sampleFG, Sample sampleBG) throws Exception {
        ModelBasedClassifier mbc = new ModelBasedClassifier(shm, myBG);
        MeasureParameters mp = new MeasureParameters(true, .99, .99, .99);
        ResultSet rs = mbc.evaluate(mp, true, sampleFG, sampleBG);
        return rs.getResults();
    }
	
    /**
     * @param shm
     * @param myBG
     * @param sampleFG
     * @param sampleBG
     * @return double array with Classification Rate at 0, Maximal Classification Rate, Area under ROC and Area under
     *         PR.
     * @throws Exception
     */
    public static double[] performClassificationTest(AbstractModel shm, AbstractModel myBG,
            Sample sampleFG,
            Sample sampleBG) throws Exception {
        ModelBasedClassifier mbc = new ModelBasedClassifier(shm, myBG);
        MeasureParameters mp = new MeasureParameters(true, .99, .99, .99);
        ResultSet rs = mbc.evaluate(mp, true, sampleFG, sampleBG);

        double dynamicCR = getDynamicClassificationRate(shm, myBG, sampleFG, sampleBG);
        System.out.println("ClassificationRate= " + (Double) (rs.getResultAt(0).getResult()));
        System.out.println("DynamicClassificationRate= " + dynamicCR);
        System.out.println("Area under ROC = " + (Double) (rs.getResultAt(7).getResult()));
        System.out.println("Area under PR = " + (Double) (rs.getResultAt(8).getResult()));
        
        return new double[] {
                (Double) rs.getResultAt(0).getResult(),
                dynamicCR,
                (Double) rs.getResultAt(7).getResult(),
                (Double) rs.getResultAt(8).getResult()
        };
    }

    /**
     * This method estimates the area under the given curve
     */
    public static double estimateArea(double[] x, double[] y) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("x and y must have the same length.");
        }

        double area = (x[0] - 0) * ((y[0] - 0) / 2 + 0);
        for (int i = 1; i < x.length; i++) {
            area += (x[i] - x[i - 1]) * ((y[i] - y[i - 1]) / 2 + y[i - 1]);
        }
        return area;
    }

    /**
     * This method computes the log-likelihood differences for the given sample and the two given models.
     * 
     * @param fg
     *            foreground model
     * @param bg
     *            background model
     * @param sample
     *            contains the sequence to work on
     * @return array of differences
     */
    public static double[] getDifferences(AbstractModel fg, AbstractModel bg, Sample sample) throws Exception {
        double[] fgScore = (fg.getLogProbFor(sample));
        double[] bgScore = (bg.getLogProbFor(sample));
        double[] diffs = new double[bgScore.length];
        for (int i = 0; i < diffs.length; i++) {
            diffs[i] = fgScore[i] - bgScore[i];
        }
        return diffs;
    }

    /**
     * Calculates the classification rate for the given threshold. All sequences whose log-likelihood difference is over
     * this threshold are assumed to be foreground, else background.
     */
    public static double getDynamicClassificationRate(
            AbstractModel fg, AbstractModel bg,
            PhyloSample sampleFG, PhyloSample sampleBG, double thresh) throws Exception {
        double[] fgDiffs = getDifferences(fg, bg, sampleFG);
        double[] bgDiffs = getDifferences(fg, bg, sampleBG);

        Arrays.sort(fgDiffs);
        Arrays.sort(bgDiffs);

        double TP = 0, FP = 0, TN = 0, FN = 0;
        for (int i = 0; i < fgDiffs.length; i++) {
            if (fgDiffs[i] > thresh) {
                TP++;
            } else {
                FN++;
            }
        }
        for (int i = 0; i < bgDiffs.length; i++) {
            if (bgDiffs[i] > thresh) {
                FP++;
            } else {
                TN++;
            }
        }
        return (TP + TN) / (TP + TN + FP + FN);
    }

    /**
     * Calculates the best classification rate. All sequences whose log-likelihood difference is over
     * this threshold are assumed to be foreground, else background
     */
    public static double getDynamicClassificationRate(
            AbstractModel fg,
            AbstractModel bg,
            Sample sampleFG,
            Sample sampleBG) throws Exception {
        double[] fgDiffs = getDifferences(fg, bg, sampleFG);
        double[] bgDiffs = getDifferences(fg, bg, sampleBG);

        Arrays.sort(fgDiffs);
        Arrays.sort(bgDiffs);

        double TP = 0, FP = 0, TN = 0, FN = 0;
        double maxRate = 0, rate = 0;
        for (int i = fgDiffs.length - 1; i >= 0; i--) {
            TP = fgDiffs.length - i;
            FN = i;

            FP = 0;
            TN = 0;
            for (int k = 0; k < bgDiffs.length; k++) {
                if (bgDiffs[k] > fgDiffs[i]) {
                    FP++;
                } else {
                    TN++;
                }
            }
            rate = (TP + TN) / (TP + TN + FP + FN);
            if (rate > maxRate) {
                maxRate = rate;
            }
        }
        return maxRate;
    }

    public static double determineThresholdForSpecifity(double specifity, AbstractModel model, Sample sample)
            throws Exception {
        int size = 0;
        for (int i = 0; i < sample.getNumberOfElements(); i++) {
            for (int u = 0; u < sample.getElementAt(i).getLength() - model.getLength() + 1; u++) {
                size++;
            }
            // size += sample.getElementAt(i).getLength() - model.getLength() + 1;
        }

        double[] scores = new double[size];
        int s = 0;
        for (int i = 0; i < sample.getNumberOfElements(); i++) {
            for (int u = 0; u < sample.getElementAt(i).getLength() - model.getLength() + 1; u++) {
                scores[s++] = model.getLogProbFor(sample.getElementAt(i), u, u + model.getLength() - 1);
            }
        }
        Arrays.sort(scores);
        return scores[(int) (scores.length * specifity)];
    }

    /**
     * This method determines for each position of each sequence in the given {@link Sample} the log probability for the
     * given {@link AbstractModel}
     * 
     * @param model
     * @param sample
     * @return an array of arrays of scores.
     * @throws Exception
     */
    public static double[][] getScores(AbstractModel model, Sample sample) throws Exception {
        if (model == null || model.getLength() == 0) {
            throw new IllegalArgumentException(
                    "Cant calculate scores for inhomogenous model");
        }

        // to calculate progress we need the size
        int sumPositions = 0;
        for (int i = 0; i < sample.getNumberOfElements(); i++) {
            sumPositions += sample.getElementAt(i).getLength() - model.getLength() + 1;
        }

        double[][] scores = new double[sample.getNumberOfElements()][];
        int count = 0; int last = 0;
        for (int i = 0; i < scores.length; i++) {
            // if sequence is shorter than motif, no scores can be calculated -> scores[i].length=0 
            scores[i] = new double[ Math.max(0, sample.getElementAt(i).getLength() - model.getLength() + 1)];
            for (int u = 0; u < sample.getElementAt(i).getLength() - model.getLength() + 1; u++) {
                scores[i][u] = model.getLogProbFor(sample.getElementAt(i), u, u + model.getLength() - 1);
                count ++;
                if(Math.ceil(count * 100. / sumPositions) > last) {
                    last = (int) Math.ceil(count * 100. / sumPositions);
                    System.out.print(last%10==0 ? "|" :".");
                }
            }
        }
        System.out.println();
        return scores;
    }

    /**
     * Determines a vector of sensitivies for the given specifities. The threshold for each specifity is calculated for
     * all positions in negativeScores. But the sensitity is determined on sequence level. I.e. a true positive is a
     * sequence with a binding site with a score larger then the determined threshold.
     * 
     * @param positiveScores
     * @param negativeScores
     * @param specifities
     * @return a vector of specifities
     */
    public static double[] getSensitivitiesPerSeqForSpecifitiesPerPos(
            double[][] positiveScores,
            double[][] negativeScores,
            double[] specifities) {

        double[] thresholds = determineThresholdsForSpecifitiesPerPos(negativeScores, specifities);
        double[] sensitivities = new double[specifities.length];
        for (int s = 0; s < thresholds.length; s++) {
            sensitivities[s] = determineSensitivityPerSeq(thresholds[s], positiveScores);
        }
        return sensitivities;
    }

    public static double[] determineThresholdsForSpecifitiesPerPos(double[][] negativeScores, double[] specifities) {
        // linearize negativeScores
        int sum = 0;
        for (int i = 0; i < negativeScores.length; i++) {
            sum += negativeScores[i].length;
        }

        double[] negScores = new double[sum];
        for (int i = 0, k = 0; i < negativeScores.length; i++) {
            for (int u = 0; u < negativeScores[i].length; u++) {
                negScores[k++] = negativeScores[i][u];
            }
        }
        Arrays.sort(negScores);

        double[] thresholds = new double[specifities.length];
        for (int s = 0; s < specifities.length; s++) {
            thresholds[s] = negScores[(int) ((negScores.length - 1) * specifities[s])];
        }
        return thresholds;
    }

    public static double determineSensitivityPerSeq(double threshold, double[][] positiveScores) {
        double truePos = 0;
        for (int i = 0; i < positiveScores.length; i++) {
            boolean isTruePos = false;
            for (int u = 0; u < positiveScores[i].length && !isTruePos; u++) {
                if (positiveScores[i][u] > threshold) {
                    truePos++;
                    isTruePos = true;
                }
            }
        }
        return truePos / positiveScores.length;
    }

    /**
     * Determines all positiones in positiveScores which have a score larger then the given threshold.
     * 
     * @param threshold
     * @param positiveScores
     * @return the set of predicted positions
     */
    public static HashSet<IntPair> determinePositionSet(double threshold, double[][] positiveScores) {
        HashSet<IntPair> set = new HashSet<IntPair>();
        for (int i = 0; i < positiveScores.length; i++) {
            boolean isTruePos = false;
            for (int u = 0; u < positiveScores[i].length && !isTruePos; u++) {
                if (positiveScores[i][u] > threshold) {
                    set.add(new IntPair(i, u));
                }
            }
        }
        return set;
    }

    public static double determineSensitivity(double threshold, AbstractModel model, Sample sample) throws Exception {
        double[][] scores = getScores(model, sample);
        return determineSensitivityPerSeq(threshold, scores);
    }

    public static int[] determineTruePositives(double threshold, AbstractModel model, Sample sample) throws Exception {
        IntList truePositives = new IntList();
        for (int i = 0; i < sample.getNumberOfElements(); i++) {
            boolean isTruePos = false;
            for (int u = 0; u < sample.getElementAt(i).getLength() - model.getLength() + 1 && !isTruePos; u++) {
                double score = model.getLogProbFor(sample.getElementAt(i), u, u + model.getLength() - 1);
                if (score > threshold) {
                    truePositives.add(i);
                    isTruePos = true;
                }
            }
        }
        return truePositives.toArray();
    }
}
