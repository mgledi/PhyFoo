package tools;

import java.io.File;
import java.io.IOException;

import javax.naming.OperationNotSupportedException;

import de.jstacs.NonParsableException;
import de.jstacs.NotTrainedException;
import de.jstacs.io.FileManager;
import de.jstacs.models.Model;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import models.PhyloPreparedAbstractModel;
import util.R;

public class PrintMotifFromSHM {
	public static void main(String[] args) throws NonParsableException, IOException, OperationNotSupportedException, NotTrainedException {
		String filename = args[0];
		String tfbs = args[1];
		StringBuffer xml = FileManager.readFile(new File(filename));
		SingleHiddenMotifMixture shm = new SingleHiddenMotifMixture(xml);

		Model model = shm.model[0];
		if (model instanceof StrandModel) {
			model = ((StrandModel) model).model[0];
		}

		System.out.println(String.format("motifs[['%s']][['condProbs%s']] = %s", tfbs, ((PhyloPreparedAbstractModel) model).getOrder(),
		        R.linearizeMatrix(((PhyloPreparedAbstractModel) model).getCondProbs())));

		/*
		 * System.out.println("Motif Model"); Util.visualizeMatrix(PWMUtil.getPWM(shm.model[0]), 3, true);
		 * 
		 * if (shm.model[0] instanceof PhyloBayesModel) { PhyloBayesModel model = ((PhyloBayesModel) shm.model[0]); System.out.println("Temperature: " +
		 * model.getBNH().getVirtualTree(0).TEMPERATURE); System.out.println("Topology: " + model.getBNH().getVirtualTree(0).getNewickString()); for (int o = 0;
		 * o < model.getBNH().getVirtualTree(0).numberOfLeafs; o++) { System.out.println("" + R.linearizeMatrix(model.getCondProbsPerSpecies(o)) +
		 * "; #R SPECIES_" + o); } }
		 * 
		 * if (shm.algorithmHasBeenRun()) { System.out.println("Best LogLikelihood: " + shm.getScoreForBestRun()); } System.out.println("ZOOPS: " +
		 * shm.getWeights()[0]);
		 * 
		 * if (shm.model[1] instanceof PhyloPreparedAbstractModel) { System.out.println("Flanking Model"); Util.visualizeMatrix(((PhyloPreparedAbstractModel)
		 * shm.model[1]).getCondProbs(), 3, true); } System.out.println("v1 = " + R.linearizeMatrix(PWMUtil.getPWM(shm.model[0])) + "; #R PWM SOLID");
		 */
	}

}
