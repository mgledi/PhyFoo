package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.utils.DoubleList;
import models.PhyloBayesModel;
import projects.dispom.PFMComparator;
import projects.dispom.PFMComparator.NormalizedEuclideanDistance;

public class EvaluationUtil {

    /** Reads the given file into a double array. Each line represents a double */
    public static double[] toVector(String filename) throws IOException {
        File f = new File(filename);
        BufferedReader br = new BufferedReader(new FileReader(f));
        String line;
        DoubleList v = new DoubleList();
        while ((line = br.readLine()) != null) {
            v.add(Double.valueOf(line));
        }
        br.close();
        return v.toArray();
    }

    /** generates a string from the given double array. The seperator can be given. */
    public static String toString(double[] d, int roundDec, String sep) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < d.length; i++) {
            if (i > 0 && i < d.length) {
                sb.append(sep);
            }
            sb.append(Util.round(d[i], roundDec));
        }
        return sb.toString();
    }

    public static double[][] toMatrix(String filename, int size) throws IOException {
        File f = new File(filename);
        BufferedReader br = new BufferedReader(new FileReader(f));
        double[][] m = new double[size][size];
        String line;
        for (int i = 0; i < size; i++) {
            line = br.readLine().trim();
            String[] vals = line.split("\t");
            for (int j = 0; j < size; j++) {
                m[i][j] = Double.valueOf(vals[j]);
            }
        }
        br.close();
        return m;
    }

    /**
     * Adds B to A and returns a pointer to A
     * 
     * @param A
     * @param B
     */
    public static void add(double[][] A, double[][] B) {
        for (int i = 0; i < B.length; i++) {
            for (int j = 0; j < B[i].length; j++) {
                A[i][j] += B[i][j];
            }
        }
    }

    /**
     * Calculates euclidean distance, Kullback Leibler Divergence and PersonCorrelationCoefficient between the two given
     * Models. Only the model at position zero is taken into account.
     * 
     * @param learnedSHM
     * @param origSHM
     * 
     * @return array of euclidean distance, Kullback Leibler Divergence and PersonCorrelationCoefficien
     */
    public static double[] compareLearnedWithOriginalPWM(SingleHiddenMotifMixture learnedSHM,
            SingleHiddenMotifMixture origSHM) {
        PhyloBayesModel origMotif = (PhyloBayesModel) origSHM.model[0];
        PhyloBayesModel learnedMotif = (PhyloBayesModel) learnedSHM.model[0];
        return compareLearnedWithOriginalPWM(learnedMotif, origMotif);
    }

    /**
     * Calculates euclidean distance, Kullback Leibler Divergence and PersonCorrelationCoefficient between the two given
     * {@link PhyloBayesModel}s and returns all three values in a double array.
     * 
     * @param learnedMotif
     * @param origMotif
     * 
     * @return array of euclidean distance, Kullback Leibler Divergence and PersonCorrelationCoefficien
     */
    public static double[] compareLearnedWithOriginalPWM(PhyloBayesModel learnedMotif, PhyloBayesModel origMotif) {
        double[][] pfm1 = origMotif.getProbabilites();
        double[][] pfm2 = learnedMotif.getProbabilites();
        return compareLearnedWithOriginalPWM(pfm1, pfm2);
    }

    /**
     * Calculates euclidean distance, Kullback Leibler Divergence and PersonCorrelationCoefficient between the two given
     * motifs and returns all three values in a double array.
     * 
     * @param pfm1
     *            first motif
     * @param pfm2
     *            second motif
     * @return array of euclidean distance, Kullback Leibler Divergence and PersonCorrelationCoefficien
     */
    public static double[] compareLearnedWithOriginalPWM(double[][] pfm1, double[][] pfm2) {
        NormalizedEuclideanDistance euclidMeasure = new NormalizedEuclideanDistance();
        PFMComparator.SymmetricKullbackLeiblerDivergence klMeasure = new PFMComparator.SymmetricKullbackLeiblerDivergence(
                0);
        PFMComparator.OneMinusPearsonCorrelationCoefficient corMeasure = new PFMComparator.OneMinusPearsonCorrelationCoefficient();
        double euclid = euclidMeasure.compare(pfm1, pfm2, pfm1.length / 2);
        double kl = klMeasure.compare(pfm1, pfm2, pfm1.length / 2);
        double cor = corMeasure.compare(pfm1, pfm2, pfm1.length / 2);
        return new double[] { euclid, kl, cor };
    }

    /**
	 * Calculates the information content of the given matrix. The returned value equals the sum of the heights of all positions within a sequence logo.
	 * 
	 * @param motif
	 * @return information content
	 */
    public static double determineIC(double[][] motif) {
        double sumH = 0;
        for (int m = 0; m < motif.length; m++) {
            sumH += determineIC(motif[m]);
        }
        return sumH;
    }
    
    public static double determineIC(double[] motifPosition) {
        double H = 2;
        for (int a = 0; a < 4; a++) {
            H += motifPosition[a] * Math.log(motifPosition[a]) / Math.log(2);
        }
        return H;
    }
}
