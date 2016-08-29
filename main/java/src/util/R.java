package util;

public class R {

    public static String getVector(int[] V) {
        StringBuffer sb = new StringBuffer();
        sb.append("c(");
        for (int a = 0; a < V.length; a++) {
            if (a != 0) {
                sb.append(", ");
            }
            sb.append(V[a]);
        }
        sb.append(")");
        return sb.toString();
    }
    public static String getVector(double[] V) {
        return getVector(V, 4);
    }
    
    public static String getVector(double[] V, int dec) {
        StringBuffer sb = new StringBuffer();
        sb.append("c(");
        for (int a = 0; a < V.length; a++) {
            if (a != 0) {
                sb.append(", ");
            }
            sb.append(Util.round(V[a], dec));
        }
        sb.append(")");
        return sb.toString();
    }
    
    public static String linearizeMatrix(double[][] m) {
        StringBuffer sb = new StringBuffer();
        sb.append("c(");
        for (int i = 0; i < m.length; i++) {
			for (int a = 0; a < m[i].length; a++) {
                if (a != 0 || i != 0) {
                    sb.append(", ");
                }
                sb.append(Util.round(m[i][a], 4));
            }
        }
        sb.append(")");
        return sb.toString();
    }
}