package tools;

import java.io.File;
import java.io.IOException;

import de.jstacs.NonParsableException;
import de.jstacs.io.FileManager;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import models.PhyloPreparedAbstractModel;
import projects.dispom.PFMComparator;
import util.PWMUtil;
import util.Util;

/**
 * This is a small tool to compare two motifs by some distance measuers. Its input should be {@link SingleHiddenMotifMixture} models, that encapsulate a
 * {@link PhyloPreparedAbstractModel}.
 * 
 * @author Martin Nettling
 * 
 */
public class CompareMotifs {
    public static void main(String[] args) throws IOException, NonParsableException {
        double[][] m1 = getMotif(args[0]);
        double[][] m2 = getMotif(args[1]);
        double[][] m2_rev = PWMUtil.getReversePWM(m2);
        
        PFMComparator.NormalizedEuclideanDistance distance = new PFMComparator.NormalizedEuclideanDistance();
        PFMComparator.SymmetricKullbackLeiblerDivergence kull = new PFMComparator.SymmetricKullbackLeiblerDivergence(0);
        
        System.out.println(Util.visualizeMatrixComparative(m1, m2, 4));
        
        
        System.out.println("Euclidean(0): " + distance.getDistance(m1, m2, 0));
        System.out.println("KL divergence(0): " + kull.getDistance(m1, m2, 0));
        
        double eucl = Math.min(distance.compare(m1, m2, m1.length / 2 + 1), distance.compare(m1, m2_rev, m1.length / 2 + 1));
        System.out.println("Euclidean(X): " + eucl);
        double kl_val = Math.min(kull.compare(m1, m2, m1.length / 2 + 1), kull.compare(m1, m2_rev, m1.length / 2 + 1));
        System.out.println("KL divergence(X): " + kl_val);
    }
    
    private static double[][] getMotif(String filename) throws IOException, NonParsableException {
        StringBuffer xml = FileManager.readFile(new File(filename));
        SingleHiddenMotifMixture shm = new SingleHiddenMotifMixture(xml);
        return PWMUtil.getPWM(shm.model[0]);
    }
}
