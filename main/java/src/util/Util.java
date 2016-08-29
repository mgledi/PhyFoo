package util;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.data.Sample;
import de.jstacs.io.FileManager;
import de.jstacs.models.AbstractModel;
import de.jstacs.models.Model;
import de.jstacs.models.discrete.inhomogeneous.BayesianNetworkModel;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.LearningType;
import de.jstacs.models.discrete.inhomogeneous.StructureLearner.ModelType;
import de.jstacs.models.discrete.inhomogeneous.parameters.BayesianNetworkModelParameterSet;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.utils.Normalisation;
import de.jstacs.utils.random.RandomNumberGenerator;
import models.PhyloBackground;
import models.PhyloBayesModel;
import models.PhyloPreparedAbstractModel;

public class Util {
    public static long randomSeed = System.currentTimeMillis();
    static Random r = new Random();

    /** gets all species from a newickstring in specific order */
    public static ArrayList<String> getOrderedArrayListFromTree(String newickString) {
        ArrayList<String> al = new ArrayList<String>();
        // clean NewickString
        newickString = newickString.replaceAll("\\s", "");
        // remove distances
        newickString = newickString.replaceAll(":[\\d]*[\\.][\\d]+", "");
        newickString = newickString.replaceAll(":[\\d]+", ":");
        // remove brackets
        newickString = newickString.replaceAll("[\\[\\]()\\\\;]", "");
        String[] ss = newickString.split(",");
        for (String s : ss) {
            al.add(s);
        }
        return al;
    }

    /** rundet ein Double auf dec Nachkommastellen */
    public static double fastRound(double val, int dec) {
        return Math.round(val * dec) / (double) dec;
    }

    /** rundet ein Double auf dec Nachkommastellen */
    public static double round(double val, int dec) {
        BigDecimal d = BigDecimal.valueOf(val);
        d = d.setScale(dec, RoundingMode.HALF_UP);
        return d.doubleValue();
    }

    /** berechnet den euklidischen Abstand zwischen zwei Arrays */
    public static double calcEuclid(double[] A, double[] B) {
        double sum = 0;
        for (int i = 0; i < A.length; i++) {
            sum += Math.pow(A[i] - B[i], 2);
        }
        return Math.sqrt(sum);
    }

    /** berechnet die Summe des �bergebenen arrays */
    public static double sum(double[] A) {
        double sum = 0;
        for (int i = 0; i < A.length; i++) {
            sum += A[i];
        }
        return sum;
    }

    /** zieht aus dem �bergebenen Array einen Index */
    public static int takeFromUniformDistr(double[] A) {
        double r = getRandomDouble(0, 1);
        double sum = 0;
        int i = 0;
        for (i = 0; i < A.length; i++) {
            if (r <= sum) {
                break;
            }
            sum += A[i];
        }
        return i - 1;
    }

    public static double[][] getRandomStochMatrix(int m, int n) {
        double[][] A = new double[m][];
        for (int i = 0; i < m; i++) {
            A[i] = getRandomStochVector(n);
        }
        return A;
    }

    /** erzeugt einen zuf�lligen Vector der Gr��e size */
    static RandomNumberGenerator rand = new RandomNumberGenerator(randomSeed);

    public static double[] getRandomStochVector(int size) {
        return getRandomStochVector(size, 1);
    }

    public static double[] getRandomStochVector(int size, double alpha) {
        double[] v = new double[size];
        for (int i = 0; i < size; i++) {
            v[i] = rand.nextGamma(alpha, 1);
        }
        Normalisation.sumNormalisation(v);
        return v;
    }

    public static double[][] getEqualStochMatrix(int m, int n) {
        double[][] A = new double[m][];
        for (int i = 0; i < m; i++) {
            A[i] = getEqualStochVector(n);
        }
        return A;
    }

    /** erzeugt einen gleichverteilten Vector der Gr��e size */
    public static double[] getEqualStochVector(int size) {
        double[] v = new double[size];
        Arrays.fill(v, 1);
        Normalisation.sumNormalisation(v);
        return v;
    }

    /** gibt eine Zufallszahl inklusive der Grenzen zur�ck */
    public static double getRandomDouble(double start, double end) {
        double t = r.nextDouble();
        t = t * (end - start) + start;
        return t;
    }

    /** gibt eine Zufallszahl inklusive der Grenzen zur�ck */
    public static int[] getRandomIntVector(int length, int start, int end) {
        int[] v = new int[length];
        for (int i = 0; i < length; i++) {
            v[i] = getRandomInt(start, end);

        }
        return v;
    }

    /** gibt eine Zufallszahl inklusive der Grenzen zur�ck */
    public static int getRandomInt(int start, int end) {
        int t = Math.abs(r.nextInt());
        t = t % (end + 1 - start) + start;
        return t;
    }

    public static double[] getGaussDistr(int vlength, double m, double s) {
        double[] norm = new double[vlength];

        for (int x = 0; x < vlength; x++) {
            norm[x] = 1 / (s * Math.sqrt(2 * Math.PI)) * Math.exp(-0.5 * Math.pow(((x - m) / s), 2));
        }
        Normalisation.sumNormalisation(norm);
        return norm;
    }

    /**
     * zieht aus einer Gaussverteilung mit der Mittelwert und Standardabweichung s
     */
    public static double getRandomGauss(double m, double s) {
        return (r.nextGaussian() * s + m);
    }

    /** kopiert das �bergebene array (1D) */
    public static double[] arraycopy(double[] a) {
        if (a == null) {
            return null;
        }

        double[] z = new double[a.length];
        z = Arrays.copyOf(a, a.length);
        return z;
    }

    /** kopiert das �bergebene array (2D) */
    public static double[][] arraycopy(double[][] a) {
        double[][] z = new double[a.length][];

        for (int i = 0; i < a.length; i++) {
            z[i] = arraycopy(a[i]);
        }
        return z;
    }

    /** kopiert das �bergebene array (3D) */
    public static double[][][] arraycopy(double[][][] a) {
        double[][][] z = new double[a.length][][];
        for (int i = 0; i < a.length; i++) {
            z[i] = arraycopy(a[i]);
        }
        return z;
    }

    /** kopiert das �bergebene array (1D) */
    public static int[] arraycopy(int[] a) {
        if (a == null) {
            return null;
        }

        int[] z = new int[a.length];
        z = Arrays.copyOf(a, a.length);
        return z;
    }

    /** kopiert das �bergebene array (2D) */
    public static int[][] arraycopy(int[][] a) {
        int[][] z = new int[a.length][];

        for (int i = 0; i < a.length; i++) {
            z[i] = arraycopy(a[i]);
        }
        return z;
    }

    /** kopiert das �bergebene array (3D) */
    public static int[][][] arraycopy(int[][][] a) {
        int[][][] z = new int[a.length][][];
        for (int i = 0; i < a.length; i++) {
            z[i] = arraycopy(a[i]);
        }
        return z;
    }

    /** visualisiert eine Matrix */
    public static String visualizeMatrix(byte[][] A, int round, boolean out) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[i].length; j++) {
                if (round == -1) {
                    sb.append(A[i][j] + "\t");
                } else {
                    sb.append(util.Util.round(A[i][j], 4) + "\t");
                }
            }
            sb.append("\n");
        }
        if (out)
            System.out.println(sb);
        return sb.toString();
    }

    /** visualisiert eine Matrix */
    public static String visualizeMatrix(int[][] A, boolean out) {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < A.length; i++) {
            if (A[i] != null) {
                for (int j = 0; j < A[i].length; j++) {
                    sb.append(A[i][j] + "\t");
                }
            } else {
                sb.append("---------");
            }
            sb.append("\n");
        }
        if (out)
            System.out.println(sb);
        return sb.toString();
    }

    /** visualisiert eine Matrix und gibt sie aus */
    public static String visualizeMatrix(double[][] A) {
        return visualizeMatrix(A, 4, true);
    }

    /** visualisiert eine Matrix */
    public static String visualizeMatrixComparative(double[][] M, double[][] N, int precision) {
        StringBuffer sb = new StringBuffer();
        if (M.length != N.length) {
            throw new IllegalArgumentException("Matrices must be of equal size. M=" + M.length + " N=" + N.length);
        }
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[i].length; j++) {
                if (precision == -1) {
                    sb.append((M[i][j]) + "\t");
                } else {
                    sb.append(round((M[i][j]), precision) + "\t");
                }
            }
            sb.append("\t||\t");
            for (int j = 0; j < N[i].length; j++) {
                if (precision == -1) {
                    sb.append((N[i][j]) + "\t");
                } else {
                    sb.append(round((N[i][j]), precision) + "\t");
                }
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    /** visualisiert eine Matrix */
    public static String visualizeMatrix(double[][] M, int precision, boolean out) {
        StringBuffer sb = new StringBuffer();

        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[i].length; j++) {
                if (precision == -1) {
                    sb.append((M[i][j]) + "\t");
                } else {
                    sb.append(round((M[i][j]), precision) + "\t");
                }
            }
            sb.append("\n");
        }
        if (out) {
            System.out.println(sb);
        }
        return sb.toString();
    }

    /** mach ein Array zu einem String */
    public static String array2string(String sep, double[] A) {
        StringBuffer tmp = new StringBuffer();
        for (int i = 0; i < A.length; i++) {
            if (i == 0) {
                tmp.append("" + String.valueOf(A[i]));
            } else {
                tmp.append(sep + "" + String.valueOf(A[i]));
            }
        }
        return tmp.toString();
    }

    /** mach ein Array zu einem String */
    public static String array2string(String sep, int[] A) {
        StringBuffer tmp = new StringBuffer();
        for (int i = 0; i < A.length; i++) {
            if (i == 0) {
                tmp.append("" + String.valueOf(A[i]));
            } else {
                tmp.append(sep + "" + String.valueOf(A[i]));
            }
        }
        return tmp.toString();
    }

    /** zieht aus einem Vector a */
    public static int getRandomIndex(double[] a, double rand) {
        for (int i = 0; i < a.length; i++) {
            rand -= a[i];

            if (rand <= 0) {
                return i;
            }
        }
        return a.length - 1;
    }

    public static double skalar(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    /**
     * Instantiates new PhyloBayes model with the given tree. All other parameters are moved from the old model.
     * 
     * @return
     * @throws Exception
     */
    public static PhyloBayesModel resetTree(PhyloBayesModel pbm, String newTree) throws Exception {
        PhyloBayesModel newPBM = new PhyloBayesModel(pbm.getLength(), pbm.getOrder(), newTree,
                pbm.getAlphabetContainer());
        for (int i = 0; i < pbm.getLength(); i++) {
            newPBM.getBNH().getVirtualTree(i).getEvolModel()
                    .setStatDistr(pbm.getBNH().getVirtualTree(i).getEvolModel().getStatDistr());
            newPBM.getBNH().getVirtualTree(i).initParametersFromGF();
        }
        return newPBM;
    }

    /**
     * Instantiates new PhyloBackground model with the given tree. All other parameters are moved from the old model.
     * 
     * @return
     * @throws Exception
     */
    public static PhyloBackground resetTree(PhyloBackground bg, String newTree) throws Exception {
        PhyloBackground newBG = new PhyloBackground(bg.getOrder(), newTree, bg.getAlphabetContainer());
        for (int i = 0; i <= bg.getOrder(); i++) {
            newBG.getBNH().getVirtualTree(i).getEvolModel()
                    .setStatDistr(bg.getBNH().getVirtualTree(i).getEvolModel().getStatDistr());
            newBG.getBNH().getVirtualTree(i).initParametersFromGF();
        }
        return newBG;
    }

    /** Returns the best model from the given directory */
    public static SingleHiddenMotifMixture getBestModel(String dir, String regex)
            throws OperationNotSupportedException, NotTrainedException, NonParsableException, IOException {
        String[] files = new File(dir).list(new RegexFileNameFilter(regex));
        double maxScore = Double.NEGATIVE_INFINITY;
        double score;
        SingleHiddenMotifMixture bestShm = null;
        for (String file : files) {
            StringBuffer xmlScore = FileManager.readFile(new File(dir + "/" + file));
            SingleHiddenMotifMixture shm = new SingleHiddenMotifMixture(xmlScore);
            score = shm.getScoreForBestRun();
            if (score > maxScore) {
                maxScore = score;
                bestShm = shm;
            }
        }
        return bestShm;
    }

    private static class RegexFileNameFilter implements FilenameFilter {
        Pattern p;

        RegexFileNameFilter(String regex) {
            super();
            p = Pattern.compile(regex);
        }

        @Override
        public boolean accept(File dir, String name) {
            Matcher m = p.matcher(name);
            if (m.matches()) {
                return true;
            }
            return false;
        }
    }

    public static double[][] getReverseFS81(double[] pi_A, double[] pi_Y, double[][] F_pi) {
        double[][] m_r = new double[pi_A.length][pi_A.length];
        for (int r = 0; r < pi_A.length; r++) {
            for (int c = 0; c < pi_A.length; c++) {
                m_r[c][r] = pi_A[r] * F_pi[r][c] / pi_Y[c];
            }
        }
        return m_r;
    }

    public static double[] times(double[] vector, double[][] matrix) {
        double[] result = new double[vector.length];
        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < matrix[0].length; col++) {
                result[col] += vector[row] * matrix[row][col];
            }
        }
        return result;
    }

    /**
     * Trys to load a {@link PhyloPreparedAbstractModel} from the given modelFile. If castToNonPhyloModel is true then the model is tried to be
     * transformed to a {@link BayesianNetworkModel}. If the motifmodel is nested in a {@link StrandModel}
     * 
     * @param modelFile
     * @param dataFG
     * @param castToNonPhyloModel
     * @return
     * @throws Exception
     */
    public static AbstractModel getMotifModel(String modelFile, Sample dataFG, boolean castToNonPhyloModel)
            throws Exception {
        // ############################## load model
        StringBuffer xml = FileManager.readFile(new File(modelFile));
        Model rawModel = new SingleHiddenMotifMixture(xml).model[0];

        if (rawModel instanceof StrandModel) {
            rawModel = ((StrandModel) rawModel).model[0];
        }

        // special case: no phylogeny is wanted and model is a BayesianNetworkModel (JStacs)
        if (rawModel instanceof BayesianNetworkModel && castToNonPhyloModel) {
            return (AbstractModel) rawModel;
        }
        
        PhyloPreparedAbstractModel model = null;
        if (rawModel instanceof PhyloPreparedAbstractModel && rawModel.getLength() > 0) {
            model = ((PhyloPreparedAbstractModel) rawModel);
        } else {
            throw new IllegalArgumentException("Model type '" + rawModel.getClass() + "'can not be handled.");
        }

        if (!castToNonPhyloModel) {
            return model;
        }

        double[][] motif, lnmotif;
        motif = model.getCondProbs();
        lnmotif = model.getLnCondProbs();

        BayesianNetworkModelParameterSet params = new BayesianNetworkModelParameterSet(
                dataFG.getAlphabetContainer(),
                model.getLength(),
                1,
                "foreground model",
                ModelType.IMM, (byte)
                0,
                LearningType.ML_OR_MAP);

        BayesianNetworkModel bnm = new BayesianNetworkModel(params);

        for (int i = 0; i < dataFG.getNumberOfElements(); i++) {
            if (dataFG.getElementAt(i).getLength() < model.getLength()) {
                System.out.println("Crap " + dataFG.getElementAt(i));
            }
        }
        bnm.train(dataFG.getInfixSample(0, model.getLength())); // must be done to initialise constraints
        for (int i = 0; i < motif.length; i++) {
            bnm.constraints[i].freq = Arrays.copyOf(motif[i], motif[i].length);
            bnm.constraints[i].lnFreq = Arrays.copyOf(lnmotif[i], lnmotif[i].length);
        }
        return bnm;
    }
    
    public static int whichMax(double[] d) {
        double m = d[0];
        int idx = 0;
        for (int i = 1; i < d.length; i++) {
            if(d[i] > m ) {
                m = d[i];
                idx = i;
            }
        }
        return idx;
    }
}
