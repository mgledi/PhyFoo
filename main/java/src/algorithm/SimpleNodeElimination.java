package algorithm;

import java.util.ArrayList;

import bayesNet.BayesNet;
import bayesNet.BayesNetNode;
import io.Alphabet;

/**
 * @author Martin Nettling TODO: Check for correctness
 */
public class SimpleNodeElimination {
    BayesNet myNet;
    int step;

    double[][] genVector;

    double[][] sumVector;

    /** help vars for fast access */
    boolean[] roots;
    int[] queueIndices;
    
    /** enth�lt die Lookup tables */
    byte[][] lookup;
    int[] pows;
    
    ArrayList<BayesNetNode> nodeQueue = new ArrayList<BayesNetNode>();
    
    public SimpleNodeElimination(BayesNet myNet) {
        nodeQueue.clear();
        queueIndices = new int[myNet.numberOfNodes];
        pows = new int[10];
        for(int i=0; i < 10; i++) {
            pows[i] = (int) Math.pow(Alphabet.size, i);
        }
        this.step = 0;
        this.myNet = myNet;

        this.lookup = new byte[pows[pows.length-1]][pows.length];
        for(int val=0; val < lookup.length; val++) {
            for(int nodePos = 0; nodePos < pows.length; nodePos ++) {
                lookup[val][nodePos] = (byte) this.getObs(val, nodePos, 0);
            }
        }
        
        this.roots = new boolean[myNet.numberOfNodes];
        for(int i=0; i < myNet.numberOfNodes; i++) {
            if(myNet.getNode(i).isRoot()) {
                this.roots[myNet.getNode(i).nodeNumber] = true;
            }
        }
        
        this.sumVector = new double[this.myNet.numberOfNodes * 2 + 1][];
        this.sumVector[0] = new double[1];
        this.sumVector[0][0] = 1;

        this.genVector = new double[this.myNet.numberOfNodes * 2 + 1][];
        this.genVector[0] = new double[1];
        this.genVector[0][0] = 1;
    }

    public double getLikelihood() {
        this.step = 0;
        for (int i = this.myNet.numberOfNodes - 1; i >= 0; i--) {
            this.outSumNode(this.myNet.getNode(i));
        }
        return this.sumVector[sumVector.length - 1][0];
    }

    /** f�gt einen Knoten der Queue hinzu, WENN dieser noch nicht enthalten ist */
    private void addNodeToQueue(BayesNetNode aktNode) {
        // TODO: node look-up-table um indexOf zu vermeiden
        if (!this.nodeQueue.contains(aktNode)) {
            this.nodeQueue.add(aktNode);
            this.queueIndices[aktNode.nodeNumber] = nodeQueue.size()-1;
        }
    }

    /** entfernt einen Knoten von der Queue */
    private void removeNodeFromQueue(BayesNetNode aktNode) {
        this.nodeQueue.remove(aktNode);
        int i=0;
        for(BayesNetNode n : nodeQueue) {
            this.queueIndices[n.nodeNumber] = i++;
        }
    }

    /** gibt die Queue aus */
    @SuppressWarnings("unused")
    private void outQueue() {
        System.out.print("QUEUE: ");
        for (int i = 0; i < this.nodeQueue.size(); i++) {
            System.out.print(this.nodeQueue.get(i).getName() + " | ");
        }
        System.out.println();
    }

    /**
     * summiert den �bergebenen knoten aus den Bereits betrachteten Knoten aus, dabei k�nnen weitere Knoten der Queue
     * hinzugef�gt werden. z.B. die Eltern des Knotens
     */
    public void outSumNode(BayesNetNode aktNode) {
        this.step++; // z�hlt mit jedem auszusummierenden Knoten hoch. Wird ben�tigt um die Reihenfolge abrufen zu
                     // k�nnen in der summiert wurde
//         System.out.println(" ************** Schritt " + step + "("+ nodeQueue.size() +") ************************ ");

        int oldLength = nodeQueue.size(); // alte Queue ist in der neuen Queue komplett enthalten
        int a; // aktuell zu betrachtendes Zeichen im aktNode, tempor�re Variable
        int newSumInd; // aktueller Index im SumVector, tempor�re Variable
        double val; // tempor�r zu merkender Wert

        // alle durch aktNode betroffenen Knoten in die Queue dazu schreiben
        if (!this.nodeQueue.contains(aktNode)) {
            this.addNodeToQueue(aktNode);
        }

        for (int p = 0; p < aktNode.numberOfParents; p++) {
            this.addNodeToQueue(aktNode.Aparents[p]);
        }
        int newLength = nodeQueue.size();
        int aktNodeInd = this.queueIndices[aktNode.nodeNumber]; // index des aktuellen Knotens in der Queue
        // ####################################### anlegen des generativen Vectors
        // ##########################################
        // ##################################################################################################################
        int lengthAllNodes = (int) pows[newLength]; // Gr��e der gesamten Potenzmenge // TODO: verallgemeinern auf m�gl. Beobachtungen
        int oldHistSize = oldLength; // merke Gr��e des alten Histograms
        int newHistSize = nodeQueue.size(); // merke Gr��e des neuen Histograms

        // genVector initialisieren
        this.genVector[this.step] = new double[lengthAllNodes];

//         System.out.println("sumVector vorher (" + (step - 1) + ") " + " ("+sumVector[step-1].length+") "+Arrays.toString(sumVector[step - 1]));
//         System.out.println("genVector vorher (" + (step) + ") " + " ("+genVector[step].length+") "+Arrays.toString(genVector[step]));

        // es k�nnen nur �nderungen am genVector auftreten, wenn neue Knoten hinzukommen, bei der Wurzel brauch nur
        // aussummiert werden
        if (!this.roots[aktNode.nodeNumber]) {
            int oldSumVecInd = 0;
            int nodeInd;

            for (int l1 = 0; l1 < lengthAllNodes; l1++) {
//                a = this.getObs(l1, aktNodeInd, newLength - 1);
                a = lookup[l1][aktNodeInd];

                // bestimme Position im CondProb-array, um auf die �bergangsWK zugreifen zu k�nnen
                int[] query = new int[aktNode.numberOfParents];
                for (int p = 0; p < aktNode.numberOfParents; p++) {
                    nodeInd = this.queueIndices[aktNode.Aparents[p].nodeNumber]; 
                    //nodeInd = this.nodeQueue.indexOf(aktNode.Aparents[p]);
//                    query[p] = this.getObs(l1, nodeInd, newLength - 1);
                    query[p] = lookup[l1][nodeInd];
                }

                // bestimme Position im alten SumVectorArray
                oldSumVecInd = 0;
                // laufe dazu �ber die bekannten Knoten im VORHERIGEN Schritt
                for (int k = 0; k < oldLength; k++) {
                    // schaue in die aktuelle Tabelle, welche realisierung der Knoten hat und multpliziere dies mit den
                    // Werten aus dem VORHERIGEN histArray, da die vorherige Position ben�tigt wird
//                    oldSumVecInd += this.getObs(l1, k, newLength - 1) *  pows[oldHistSize - 1 - k];
                    oldSumVecInd += lookup[l1][k] *  pows[oldHistSize - 1 - k];
                }
                // System.out.println(aktNode.CPF.get(query, a) + " x " + this.sumVector[step-1][oldSumVecInd]);

                // Berechne die neuen Werte im Array
                this.genVector[this.step][l1] = aktNode.CPF.get(query, a) * this.sumVector[step-1][oldSumVecInd];
            }
        }
        // wurde eine Wurzel erreicht, reicht es den alten SumVector als generativen Vector zu verwenden
        else {
            this.genVector[this.step] = this.sumVector[this.step - 1];
        }

        // ####################################### aussummieren der Werte aus dem generativen Vector in den SumVector
        // #######
        // ##################################################################################################################
        // initalisiere den neuen sumVector mit der Gr��e der Potenzmenge aller Knoten dividiert mit der Alphabetgr��e,
        // da genau so viele Werte durch
        // das aussummieren entfallen werden
        this.sumVector[this.step] = new double[lengthAllNodes / Alphabet.size]; // TODO: verallgemeinern auf m�gl. Beobachtungen

        // laufe �ber alle m�glichen Werte im GenVector und summiere entsprechende miteinander
        for (int l1 = 0; l1 < lengthAllNodes; l1++) {
//            a = this.getObs(l1, aktNodeInd, newLength - 1);
            a = lookup[l1][aktNodeInd];

            newSumInd = 0;
            for (int k = 0; k < newLength; k++) {
                // da aktnode gleich aus der queue entfernt wird, muss sich das Hist ohne aktNode vorgestellt werden.
                // das hei�t aktnode wird nicht mitbetrachtet und alle Knoten nach aktNode rutschen eine Stelle vor
                if (k > aktNodeInd) {
//                    newSumInd += this.getObs(l1, k, newLength - 1) * pows[newHistSize - 1 - k] ;//this.hist[step][newHistSize - 1 - k];
                    newSumInd += lookup[l1][k] * pows[newHistSize - 1 - k] ;//this.hist[step][newHistSize - 1 - k];
                } else if (k < aktNodeInd) {
//                    newSumInd += this.getObs(l1, k, newLength - 1) * pows[newHistSize - 1 - k - 1]  ;//this.hist[step][newHistSize - 1 - k - 1];
                    newSumInd += lookup[l1][k] * pows[newHistSize - 1 - k - 1]  ;//this.hist[step][newHistSize - 1 - k - 1];
                }
            }
            // wenn AktNode beobachtet ist oder eine Wurzel, dann ist Prob gesetzt und es findet eine Multiplikation statt
            if (aktNode.props._isObserved && roots[aktNode.nodeNumber]) {
                val = this.genVector[step][l1] * aktNode.props._observation[a] * aktNode.CPF.get(0, a);
            } else if (aktNode.props._isObserved) {
                val = this.genVector[step][l1] * aktNode.props._observation[a]; //
            } else if (roots[aktNode.nodeNumber]) {
                // TODO: Wurzel kann auch beobachtet sein
                val = this.genVector[step][l1] * aktNode.CPF.get(0, a);;
            }
            // ansonsten wird ohne Multiplizieren aufsummiert
            else {
                val = this.genVector[step][l1];
            }
            this.sumVector[step][newSumInd] += val;
        }

        // markieren, dass der gerade angefasste Knoten nun aussummiert wird
        this.removeNodeFromQueue(aktNode);
        this.step++;
        // this.generateLookUpTable();
        this.sumVector[step] = this.sumVector[step - 1];
        this.genVector[step] = this.genVector[step - 1];

        this.genVector[step - 1] = null;
        this.sumVector[step - 1] = null;
//         System.out.println("sumVector(" + (this.step) + ") " + Arrays.toString(this.sumVector[this.step]));
//         System.out.println("size: " + nodeQueue.size());
//         System.out.println("genVector(" + (this.step) + ") " + Arrays.toString(this.genVector[this.step]));
//         this.outQueue();

    }

    /** 
	 */
    public int getObs(int val, int nodePos, int maxNodePos) {
        return (int) ((val / pows[nodePos]) % Alphabet.size);
//        return (int) ((val / Math.pow(Alphabet.size, nodePos)) % Alphabet.size);
//        return (int) ((val / Math.pow(Alphabet.size, maxNodePos - nodePos)) % Alphabet.size);
    }
    
    
    
    public static SimpleNodeElimination INSTANCE;
    public static void instantiate(BayesNet net) {
        if( INSTANCE  == null ) {
            INSTANCE = new SimpleNodeElimination(net);
        }
    }
}