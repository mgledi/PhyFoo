package models;

import de.jstacs.NonParsableException;
import de.jstacs.data.AlphabetContainer;
import de.jstacs.utils.Normalisation;

public abstract class AbstractBoltzmanModel extends AbstractAlignmentBasedModel {

    double[][] temperatures;
    double[][][] bendedCondProbs;
    int numberOfSpecies = 1;

    public AbstractBoltzmanModel(AlphabetContainer alphabets, byte order, int length, int species) {
        super(alphabets, order, length);
        this.numberOfSpecies = species;
    }

    public AbstractBoltzmanModel(StringBuffer xml) throws NonParsableException {
        super(xml);
    }

    public void reinitBendedCondProbs() {
        for (int o = 0; o < temperatures.length; o++) {
            for (int m = 0; m < condProb.length; m++) {
                double[] bended = new double[condProb[m].length];
                for (int a = 0; a < bended.length; a++) {
                    bended[a] = Math.pow(condProb[m][a], temperatures[o][temperatures[o].length == 1 ? 0 : m]);
                }
                Normalisation.sumNormalisation(bended); //TODO: this only works for PWMs
                bendedCondProbs[o][m] = bended;
            }
        }
    }

    public int getNumberOfTemperaturesToOptimize() {
        return temperatures.length * temperatures[0].length;
    }
    
    public double[] getLinearizedTemperatures() {
        int k=0; 
        double[] temp = new double[getNumberOfTemperaturesToOptimize()]; 
        for (int o = 0; o < temperatures.length; o++) {
            for (int m = 0; m < temperatures[o].length; m++) {
                temp[k++] = temperatures[o][m];
            }
        }
        return temp;
    }
    
    public void setLinearizedTemperatures(double[] temp) {
        int k=0; 
        for (int o = 0; o < temperatures.length; o++) {
            for (int m = 0; m < temperatures[o].length; m++) {
                temperatures[o][m] = temp[k++];
            }
        }
    }
    
    /**
     * Sets the temperature for a given species and a given position. If the temperature should be learned globally
     * position must be set to 0.
     * 
     * @param species
     * @param position
     * @param temperature
     */
    public void setTemperature(int species, int position, double temperature) {
        this.temperatures[species][position] = temperature;
    }

    /**
     * @param species
     * @param position
     * @return the temperature for the given species and the given position
     */
    public double getTemperature(int species, int position) {
        return this.temperatures[species][position];
    }

    /**
     * Accesses the number of species
     * 
     * @return the number of species
     */
    public int getNumberOfSpecies() {
        return numberOfSpecies;
    }
}
