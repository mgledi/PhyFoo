package algorithm;

import java.util.Arrays;

import bayesNet.BayesNet;
import bayesNet.BayesNetHandler;
import bayesNet.BayesNetNode;
import io.Alphabet;

/**
 * Algorithmus zur effektiven Berechnung der Likelihood eines Bayesnetzes mit beliebigen Beobachtungen
 * 
 * @author Chaos
 */
public class NodeSumming {
    private BayesNet _myNet;
    private BayesNetHandler bnh;

    /**
     * bildet die Knotennummer auf eine Algorithmus spezifische Id ab, da intern nur unbeobachtet Knoten interessant
     * sind
     */
    private int[] _nodeNumber2hiddenNodeNumber;

    /** zahl der unbeobachteten Zust�nde */
    private int _numberHiddenNodes = 0;

    private int[][] observations;

    public NodeSumming(BayesNetHandler bnh) {
        this.bnh = bnh;
        this._myNet = bnh.getNet();
        this.generateHiddenNodeLookup();
        this.generateAllObservations();
    }

    /**
     * generiert eine Abbildung der Nodenumber auf die HiddenNodenumber, so muss sp�ter nicht �ber alle sondern wirklich
     * nur �ber die unbeobachteten Knoten gelaufen werden
     */
    private void generateHiddenNodeLookup() {
        _nodeNumber2hiddenNodeNumber = new int[_myNet.numberOfNodes];

        int j = 0;
        for (int i = 0; i < _myNet.numberOfNodes; i++) {
            if (!_myNet.getNode(i).props.isObserved()) {
                _nodeNumber2hiddenNodeNumber[_myNet.getNode(i).nodeNumber] = j;
                _numberHiddenNodes++;
                j++;
            }
        }
    }

    private void generateAllObservations() {
        int size = (int) Math.pow(Alphabet.size, _numberHiddenNodes);
        
        observations = new int[size][];
        for (int i = 0; i < size; i++) {
            // hole die beobachtung zu i
            observations[i] = getLookUpVectorFromVal(i);
        }
    }

    /**
     * Calculates the Likelihood up to the given Position
     * 
     * @param end
     * @return likelihood
     */
    public double getLikelihood(int end) {
        double sum = 0;
        int size = (int) Math.pow(Alphabet.size, (bnh.getVirtualTree(0).numberOfNodes-bnh.getVirtualTree(0).numberOfLeafs) * (end+1));
        for (int i = 0; i < size; i++) {
            // hole die beobachtung zu i
//            System.out.println(Arrays.toString(observations[i]));
            sum += this.getLikelihoodForCompleteObs(observations[i], end);
        }
        return sum;
    }

    /**
     * berechnet die Likelihood zu einer vollst�ndigen Beobachtung, ACHTUNG: obs enth�lt nur die Beobachtung zu den
     * unbekannten Knoten, alle anderen Beobachtungen stehen in der Netzstruktur
     * 
     * @param obs
     * @param end
     * 
     * @return likelihood
     */
    public double getLikelihoodForCompleteObs(int[] obs, int end) {
        BayesNetNode aktNode;
        double product = 1;
        double factor;
        int[] obsParent;
        
        for (int k=0; k <= end; k++) {
            for(int i =0; i < bnh.getVirtualTree(k).numberOfNodes; i++) {
//                System.out.println("Tree " + k + " node " + i);
                aktNode = bnh.getVirtualTree(k).getNode(i);
                // System.out.println("======> " + aktNode.getName());
                obsParent = new int[aktNode.numberOfParents];
                for (int p = 0; p < aktNode.numberOfParents; p++) {
                    if (aktNode.getParent(p).props.isObserved()) {
                        obsParent[p] = aktNode.getParent(p).props.getObservedIndex();
                    } else {
                        obsParent[p] = obs[_nodeNumber2hiddenNodeNumber[aktNode.getParent(p).nodeNumber]];
                    }
                }
                // System.out.println("parent: " + Util.array2string(",",obsParent));
                
                if (aktNode.props.isObserved()) {
                    factor = aktNode.CPF.get(obsParent, aktNode.props.getObservedIndex());
                } else {
                    factor = aktNode.CPF.get(obsParent, obs[_nodeNumber2hiddenNodeNumber[aktNode.nodeNumber]]);
                }
                // System.out.println(product + " * " + factor);
                product *= factor;
            }
        }
//        for (int i = 0; i < _myNet.numberOfNodes; i++) {
//        }
        return product;
    }

    private int[] getLookUpVectorFromVal(int val) {
        // initialisiere Anfrage
        int[] query = new int[_numberHiddenNodes];
        Arrays.fill(query, -1);
        // laufe �ber all Eltern um aus dem �bergebenen Index die beobachtung zu berechnung
        for (int i = _numberHiddenNodes - 1; i >= 0; i--) {

            int size = (int) Math.pow(Alphabet.size, i);
            if (i == 0) {
                query[0] = (int) val % Alphabet.size;
            } else {
                query[i] = (int) val / size % Alphabet.size;
            }
        }
        return query;
    }
}