package algorithm;

import org.apache.log4j.Logger;

import de.jstacs.NonParsableException;
import de.jstacs.data.Sample;
import de.jstacs.models.mixture.StrandModel;
import de.jstacs.models.mixture.motif.SingleHiddenMotifMixture;
import io.FileUtil;
import models.PhyloPreparedAbstractModel;

/**
 * This class handles several iterations of an EM algorithm containing only one model. If a new optimum was found the
 * last written file is overwritten with the current model. This handler was implemented, because the functionality of
 * the internal repeat mechanism was not sufficient.
 * 
 * @author Martin Nettling
 */
public class EMHandler {

	private static Logger LOGGER = Logger.getLogger(EMHandler.class);

    String initialModel;
    private int repeats;
    private Sample data;

    private double[] finalScores;
    private String[] finalModels;
    private int iteration = 0;

    private EMHandler(SingleHiddenMotifMixture initialModel, int repeats, Sample data) {
        if (repeats < 1) {
            repeats = 1;
        }
        this.initialModel = initialModel.toXML().toString();
        this.repeats = repeats;
        this.data = data;
        finalScores = new double[repeats];
        finalModels = new String[repeats];
    }

    /**
     * This method internally instantiates an EMHandler. The given {@link SingleHiddenMotifMixture} will be trained
     * several times on the given data. If the given saveFileName is not null, after each iteration the best model is
     * stored in this file.
     * 
     * @param initialModel
     * @param repeats
     * @param data
     * @param saveFileName
     * @return this {@link EMHandler}
     * @throws Exception
     */
    public static EMHandler train(SingleHiddenMotifMixture initialModel, int repeats, Sample data, String saveFileName)
            throws Exception {
        PhyloPreparedAbstractModel motifModel = null;
        if(initialModel.model[0] instanceof PhyloPreparedAbstractModel) {
            // do nothing
        } else if (initialModel.model[0] instanceof StrandModel && ((StrandModel)initialModel.model[0]).model[0] instanceof PhyloPreparedAbstractModel) {
            // do nothing
        } else {
            if(initialModel.model[0] instanceof StrandModel) {
                throw new IllegalArgumentException("Can not handle " + ((StrandModel)initialModel.model[0]).model[0].getClass());
            } else {
                throw new IllegalArgumentException("Can not handle " + initialModel.model[0].getClass());
            }
        }
        
        EMHandler em = new EMHandler(initialModel, repeats, data);
        double best=Double.NEGATIVE_INFINITY;
        for (; em.iteration < em.repeats; em.iteration++) {
			LOGGER.info("EMHandler: Start " + (em.iteration + 1) + " of " + repeats);
            SingleHiddenMotifMixture shm = new SingleHiddenMotifMixture(new StringBuffer(em.initialModel));
            if(shm.model[0] instanceof PhyloPreparedAbstractModel) {
                motifModel =  (PhyloPreparedAbstractModel)shm.model[0];           
            } else if (shm.model[0] instanceof StrandModel && ((StrandModel)shm.model[0]).model[0] instanceof PhyloPreparedAbstractModel) {
                motifModel =  (PhyloPreparedAbstractModel)((StrandModel)shm.model[0]).model[0];           
            }
			if (repeats > 1) {
				motifModel.fillLnCondProbRandomly();
			}
            shm.train(em.data);
            em.finalModels[em.iteration] = shm.toXML().toString();
            em.finalScores[em.iteration] = shm.getScoreForBestRun();
            if(em.iteration == 0 || em.finalScores[em.iteration] > best) {
                best = em.finalScores[em.iteration];
				LOGGER.info("New best run (" + em.iteration + ") = " + em.finalScores[em.iteration]);
            }
            if (saveFileName != null) {
				FileUtil.writeFile(saveFileName, em.getBestModel().toXML().toString());
            }
        }
        return em;
    }
    
    /**
     * Returns the best model concerning the likelihood in {@link SingleHiddenMotifMixture}.
     * 
     * @return the best {@link SingleHiddenMotifMixture} learned up to now.
     * @throws NonParsableException
     */
    public SingleHiddenMotifMixture getBestModel() throws NonParsableException {
        String bestModel = finalModels[0];
        double bestScore = finalScores[0];
        for (int i = 1; i < iteration; i++) {
            if (finalScores[i] > bestScore && finalModels[i] != null) {
                bestScore = finalScores[i];
                bestModel = finalModels[i];
            }
        }
        return new SingleHiddenMotifMixture(new StringBuffer(bestModel));
    }
}
