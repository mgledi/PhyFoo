package tools;

import java.util.Arrays;
import java.util.Properties;

import classification.ClassificationUtil;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import de.jstacs.results.Result;
import de.jstacs.utils.DoubleList;
import io.FileUtil;
import io.SampleUtil;
import models.ModelUtil;
import models.PhyloBackground;
import models.PhyloBayesModel;
import models.PhyloPreparedAbstractModel;
import util.Config;

/**
 * Example configuration:
 * --basedir="e:/Dropbox/promotion/edgeLengthTest_FS81alpha/biological_tests/CTCF_gapless" <br>
 * --input.dataset.fg=%basedir%/dataset/data_1.0.fg <br>
 * --input.dataset.bg=%basedir%/dataset/data_1.0.bg <br>
 * --model.newick=(SPECIES_0:1.0,SPECIES_1:1.0,SPECIES_2:1.0,SPECIES_3:1.0,SPECIES_4:1.0) <br>
 * --input.splitseed=0 <br>
 * --input.dataset.fg.train=200 <br>
 * --input.dataset.bg.train=500 <br>
 * --input.dataset.bg.test= <br>
 * --output.results=%basedir%/dataset%input.splitseed%/results.dat
 * </code>
 */
public class Classification {
	public static void main(String[] args) throws Exception {
		Properties props = new Properties();
		Config.parseProperties(props, args, false);
		Config.replaceInternalPlaceholder(props);

		AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

		// get trained SHM to predict positions
		SingleHiddenMotifMixture bestSHM = ModelUtil.getTrainedSHM(props);
		PhyloPreparedAbstractModel bestBG = ModelUtil.getTrainedBackground(props);
		
		String resultFile = Config.getProperty(props, "output.results", false).asString();

		// get datasets
		Sample[] samples = SampleUtil.getDataSets(props, alphabet);
		Sample sampleTestFG = samples[1], sampleTrainBG = samples[2], sampleTestBG = samples[3];

		// Hack: In some cases the background was not trained correctly. reinit it, and post train
		// ((PhyloBackground)bestBG).reinitEdgelengths(((PhyloBackground)bestSHM.model[1]).getBNH().getVirtualTree(0).getNewickString());
		// ModelUtil.trainModel(bestBG, sampleTrainBG, 1e-2, 2, null);
	
		PhyloBackground.ENABLE_CACHING = false;
        PhyloBayesModel.ENABLE_CACHING = false;
        PhyloBayesModel.CLONING_ALLOWED = false;

        System.out.println("Start tests on " + sampleTestFG.getNumberOfElements() + " FG-seqs and " + sampleTestBG.getNumberOfElements() + " BG-seqs.");
//        double dynamicClassificationRate = ClassificationUtil.getDynamicClassificationRate(bestSHM, bestBG, sampleTestFG, sampleTestBG);
//        System.out.println(String.format("%s : %f", "Dynamic Classification rate", dynamicClassificationRate));
        Result[] results = ClassificationUtil.performClassification(bestSHM, bestBG, sampleTestFG, sampleTestBG);
        double classificationRate = (Double) results[0].getResult();
        System.out.println(String.format("%s : %f", results[0].getName(), classificationRate));
        double areaUnderROC = (Double) results[7].getResult();
        System.out.println(String.format("%s : %f", results[7].getName(), areaUnderROC));
        double areaUnderPR = (Double) results[8].getResult();
        System.out.println(String.format("%s : %f", results[8].getName(), areaUnderPR));
        double[][] rocCurve = (double[][])results[11].getResult();
        System.out.println(String.format("%s : %f ...", results[11].getName(), rocCurve[0][0]));
        double[][] prCurve = (double[][])results[12].getResult();
        System.out.println(String.format("%s : %f ...", results[12].getName(), prCurve[0][0]));

        StringBuilder sb = new StringBuilder();
        sb.append(classificationRate).append("\n");
//        sb.append(dynamicClassificationRate).append("\n");
        sb.append(areaUnderROC).append("\n");
        sb.append(areaUnderPR).append("\n");
		sb.append("ROC" + "[['x']]").append(Arrays.toString(getColumn(rocCurve, 0))).append("\n");
		sb.append("ROC" + "[['y']]").append(Arrays.toString(getColumn(rocCurve, 1))).append("\n");
		sb.append("PR" + "[['x']]").append(Arrays.toString(getColumn(prCurve, 0))).append("\n");
		sb.append("PR" + "[['y']]").append(Arrays.toString(getColumn(prCurve, 1))).append("\n");
        // store model and results
        System.out.println(sb);

		System.out.println("Writing data to " + resultFile);
		FileUtil.writeFile(resultFile, sb.toString());
	}
	
	private static double[] getColumn(double[][] m, int idx) {
		DoubleList column = new DoubleList();
		for (int i = 0; i < m.length; i++) {
			column.add(m[i][idx]);
		}
		return column.toArray();
	}
}
