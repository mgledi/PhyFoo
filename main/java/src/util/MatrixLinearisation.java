package util;

import de.jstacs.utils.Normalisation;

public class MatrixLinearisation {

    /**
     * wandelt einen Vector vom lambda-raum in den pi-Raum um
     * 
     * @param lambda
     * @param pi the target
     * @param alphabetSize 
     * 
     */
    public static void lambda2pi(double[] lambda, double[] pi, int alphabetSize) {
        for (int i = 0; i < lambda.length; i += alphabetSize) {
            Normalisation.logSumNormalisation(lambda, i, i + alphabetSize, pi, i);
        }
    }

    /**
     * wandelt einen Vector vom lambda-raum in den pi-Raum um
     * 
     * @param double[] lambda
     * 
     * @return double[] pi
     */
    public static double[] lambda2pi(double[] lambda, int alphabetSize) {
        double[] vec = new double[lambda.length];
        System.arraycopy(lambda, 0, vec, 0, lambda.length);
        for (int i = 0; i < lambda.length; i += alphabetSize) {
            Normalisation.logSumNormalisation(lambda, i, i + alphabetSize, vec, i);
        }
        return vec;
    }

    /**
     * wandelt eine Matrix vom lambda-raum in den pi-Raum um ACHTUNG:
     * Matrix-zeilen sollten gleich lang sein
     * 
     * @param double[][] lambda
     * 
     * @return double[][] pi
     */
    public static double[][] lambda2pi(double[][] lambda) {
        double[][] vec = new double[lambda.length][];
        for (int i = 0; i < lambda.length; i++) {
            vec[i] = lambda2pi(lambda[i], lambda[i].length);
        }
        return vec;
    }

    /**
     * wandelt einen Vector vom pi- in den Lambda-Raum um ACHTUNG: Matrix-zeilen
     * sollten gleich lang sein
     * 
     * @param double[][] pi
     * 
     * @return double[][] lambda
     */
    public static double[][] pi2lambda(double[][] pi) {
        double[][] vec = new double[pi.length][];
        for (int i = 0; i < vec.length; i++) {
            vec[i] = pi2lambda(pi[i]);
        }
        return vec;
    }

    /**
     * wandelt einen Vector vom pi- in den Lambda-Raum um
     * 
     * @param pi
     * 
     * @return double[] lambda
     */
    public static double[] pi2lambda(double[] pi) {
        double[] vec = new double[pi.length];
        for (int i = 0; i < vec.length; i++) {
            vec[i] += Math.log(pi[i]);
        }
        return vec;
    }

    /**
     * linearisiert eine Matrix: ACHTUNG: alle Zeilen sollten gleich lang sein
     * 
     * @param m zu linearisierende Matrix
     */
    public static double[] linearize(double[][] m) {
        int size = m.length * m[0].length;
        double[] l = new double[size];
        // laufe über alle Spalten
        int k = 0;
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[i].length; j++) {
                l[k] = m[i][j];
                k++;
            }
        }
        return l;
    }
    
    public static void fillMatrix(double[] v, double[][] m) {
        int k = 0;
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[i].length; j++) {
                m[i][j] = v[k++];
            }
        }
    }
    
    public static double[][] toMatrix(double[] v, int dimX, int dimY) {
        double[][] m = new double[dimY][dimX];

        int k = 0;
        for (int i = 0; i < m.length; i++) {
            for (int j = 0; j < m[i].length; j++) {
                m[i][j] = v[k++];
            }
        }
        return m;
    }
}
