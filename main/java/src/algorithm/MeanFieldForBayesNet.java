package algorithm;

import java.util.Arrays;

import bayesNet.BayesNet;
import bayesNet.BayesNetHandler;
import bayesNet.BayesNetNode;
import de.jstacs.data.Sequence;
import de.jstacs.data.Sequence.SubSequence;
import de.jstacs.data.sequences.MultiDimensionalDiscreteSequence;
import io.Alphabet;

/**
 * Calculates and optimizes the meanfield parameters of a {@link BayesNet}
 * 
 * @author Martin Nettling
 */
public class MeanFieldForBayesNet {

    /** Use this parameter for DEBUG purposes only */
    public static boolean CALC_LIKELIHOOD = false;

    /** Pointer to the BayesNetHanlder to use */
    private BayesNetHandler bnh;

    /** approximations of distributions in unobserved nodes */
    private double[] q;

    /** pointer to the underlying sequence */
    public Sequence<int[]> seq;

    /**
     * pointer to the underlying MultiDimensionalDiscreteSequence of {@link #seq} underlying sequence. If seq is a
     * {@link MultiDimensionalDiscreteSequence} than seq == parent
     */
    public MultiDimensionalDiscreteSequence parent;

    private int[] node2hiddenNodeMap;
    private int hiddenNodes;

    public MeanFieldForBayesNet(BayesNetHandler bnh, Sequence<int[]> seq) {
        this.bnh = bnh;
        this.seq = seq;
        // determine the underlying MultiDimensionalDiscreteSequence for faster lookups
        if (seq instanceof MultiDimensionalDiscreteSequence) {
            parent = (MultiDimensionalDiscreteSequence) seq;
        } else {
            parent = (MultiDimensionalDiscreteSequence) ((SubSequence<?>) seq).getParent();
        }
    }

    public void init() {
        // estimate number of unkown nodes
        hiddenNodes = 0;
        node2hiddenNodeMap = new int[bnh.getNet().numberOfNodes];
        for (int i = 0; i < this.bnh.getNet().numberOfNodes; i++) {
            if (!bnh.getNet().getNode(i).props._isObserved) {
                node2hiddenNodeMap[i] = hiddenNodes * 4;
                hiddenNodes++;
            } else {
                node2hiddenNodeMap[i] = -1;
            }
        }
        if (hiddenNodes == bnh.getNet().numberOfNodes) {
            System.out.println("+++++ Warning +++++ MeanFieldForBayesNet: There are no observed nodes.");
        }
        this.q = new double[hiddenNodes * Alphabet.size];
        Arrays.fill(this.q, 1. / Alphabet.size); // initializing uniform distributed
    }

    public void setSequence(Sequence<int[]> seq) {
        this.seq = seq;
    }

    /**
     * initialisiert die Beobachtung f�r das Netz, sollte immer vor irgendeiner Berechnung passieren
     */
    public void initObservation() {
        int offset = 0;
        if (!(seq instanceof MultiDimensionalDiscreteSequence)) {
            offset = ((SubSequence<?>) seq).start;
        }

        for (int m = 0; m < Math.min(bnh.motifLength, seq.getLength()); m++) {
            // _bnh._myTreesArray[m].setObservation(getObservation(m));
            int[] tmp = parent.containerForPhyloBayes[m + offset];
            bnh._myTreesArray[m].setObservation(tmp);
        }
    }

    /** optimizes the free energy concerning q (the meanfields) */
    public void optimizeByNormalisation() {
        if (CALC_LIKELIHOOD) {
            return;
        }
        BayesNetNode aktNode, aktChild;
        double prod = 1;
        int[] query;
        double[] newQ = new double[q.length];

        // variables for convergence
        int maxSteps = 100, minsteps = 2;
        boolean converged = false;
        int k = 0;
        double oldScore = calcFreeEnergy(), newScore;

        while (k < maxSteps && !converged || k < minsteps) { // minsteps is needed
            NODES: for (int i = 0; i < bnh.getNet().numberOfNodes; i++) { // run over all nodes
                aktNode = bnh.getNet().getNode(i); // remember actual node
                if (aktNode.props._isObserved) { // if the current node is observed, nothing must be done
                    continue NODES;
                }
                for (int a = 0; a < Alphabet.size; a++) { // run over all possible obeservations
                    // System.out.println("-------------- Behandle " + aktNode.getName() + " symb " + a);
                    newQ[node2hiddenNodeMap[i] + a] = 0;
                    // ####################### berechne Grundsumme des Knotens
                    for (int l = 0; l < aktNode.CPF._size; l++) {
                        query = aktNode.CPF._lookup[l]; // getQueryByIndex
                        prod = 1;
                        for (int p = 0; p < aktNode.numberOfParents; p++) {
                            if (!aktNode.Aparents[p].props._isObserved) { // nur das Q eines eltern ansprechen, wenn
                                                                          // dieser nicht beobachtet ist (sonst gibt es
                                                                          // kein q)
                                prod *= this.q[node2hiddenNodeMap[aktNode.Aparents[p].nodeNumber] + query[p]];
                                // System.out.println(" * " +
                                // this.q[_nodeNumber2hiddenNodeNumber[aktNode.Aparents[p].nodeNumber]][query[p]]);
                            }
                        }
                        // TODO: nur summieren wenn CPF stimmt (�berpr�fen der beobachteten eltern
                        // System.out.println(" + " + prod + " * " + Math.log(aktNode.CPF.get(l,a)));
                        newQ[node2hiddenNodeMap[i] + a] += prod * Math.log(aktNode.CPF.get(l, a));
                        // if(aktNode.CPF.get(l,a) == 0) {
                        // System.err.println("schei�e");
                        // }
                    }
                    // System.out.println(aktNode.getName() + " " + sum);

                    // Kindersumme des Knotens
                    for (int c = 0; c < aktNode.numberOfChilds; c++) {
                        aktChild = aktNode.Achildren[c];
                        // System.out.println("      ------- " + aktChild.getName());
                        for (int ac = 0; ac < Alphabet.size; ac++) // laufe �ber jede m�gliche Realisierung des
                                                                   // Kindes
                        {
                            // wenn das aktuelle kind beobachtet ist, dann nur aufsummieren, wenn die Beobachtung
                            // getroffen wurde
                            if (!aktChild.props._isObserved || aktChild.props._isObserved
                                    && aktChild.props._observedIndex == ac) {
                                for (int l = 0; l < aktChild.CPF._size; l++) {
                                    boolean trueObs = true;
                                    query = aktChild.CPF._lookup[l]; // getQueryByIndex(l);
                                    if (aktChild.props._isObserved) {
                                        prod = 1;
                                    } else {
                                        // System.out.println(Arrays.toString(_nodeNumber2hiddenNodeNumber));
                                        // System.out.println("step " + k + "  node " + aktChild.nodeNumber + " parent "
                                        // + aktNode.nodeNumber + " " +ac + " =?= " + aktChild.props._observedIndex);
                                        prod = this.q[node2hiddenNodeMap[aktChild.nodeNumber] + ac];
                                    }
                                    // laufe �ber alle Eltern des Kindes
                                    for (int p = 0; p < aktChild.numberOfParents; p++) {
                                        // �berpr�fe ob in Query die wahre Beobachtung steht
                                        if (aktChild.Aparents[p].props._isObserved
                                                && aktChild.Aparents[p].props._observedIndex != query[p]) {
                                            trueObs = false;
                                        }
                                        // nur q aufmultiplizieren, wenn es nicht gerade abgeleitet wurde, also das q
                                        // des aktuellen Knotens ist
                                        // oder der Vater beobachtet wurde (dann gibt es kein q -> nullpointer)
                                        if (aktChild.Aparents[p].nodeNumber != aktNode.nodeNumber &&
                                                !aktChild.Aparents[p].props._isObserved) {
                                            // System.out.println("      * " +
                                            // this.q[_nodeNumber2hiddenNodeNumber[aktChild.Aparents[p].nodeNumber]][query[p]]);
                                            prod *= this.q[node2hiddenNodeMap[aktChild.Aparents[p].nodeNumber]
                                                    + query[p]];
                                        }
                                    }
                                    // nur aufsummieren, wenn die Beobachtung die des elternknotens aktnode ist
                                    // TODO: nur aufsummieren wenn query = wirklichen beobachtungen der eltern
                                    // etnspricht
                                    if (query[aktChild.getParentIndex(aktNode)] == a && trueObs) {
                                        newQ[node2hiddenNodeMap[i] + a] += prod * Math.log(aktChild.CPF.get(l, ac));
                                        // if(aktChild.CPF.get(l,ac) == 0) {
                                        // System.err.println(aktChild.CPF.get(l,ac));
                                        // }
                                        // System.out.println("      " + prod + " * " + aktChild.CPF.get(l,ac));
                                    }
                                }
                            }
                        }
                    }
                    newQ[node2hiddenNodeMap[i] + a] = Math.exp(newQ[node2hiddenNodeMap[i] + a]);
                }

                double sum = 0;
                for (int a = 0; a < 4; a++) {
                    sum += newQ[node2hiddenNodeMap[i] + a];
                }
                for (int a = 0; a < 4; a++) {
                    this.q[node2hiddenNodeMap[i] + a] = newQ[node2hiddenNodeMap[i] + a] / sum;
                }
            }

            newScore = calcFreeEnergy();
            if (Math.abs(oldScore - newScore) < 1e-5) {
                converged = true;
            }
            oldScore = newScore;
            k++;
        }
    }

    /** berechnet aus dem Klasse-q und einem Virtuellen Baum die Freie Energie */
    public double calcFreeEnergy() {
        return calcFreeEnergy(0, seq.getLength() - 1);
    }

    /** berechnet aus dem Klasse-q und einem Virtuellen Baum die Freie Energie */
    NodeSumming ns;

    public double calcLogLikelihood(int start, int end) {
        // NodeSumming ns = new NodeSumming(_bnh); // TODO: buggy
        // return Math.log(ns.getLikelihood(end)) - Math.log(ns.getLikelihood(start));
        SimpleNodeElimination.instantiate(this.bnh.getNet());
        return Math.log(SimpleNodeElimination.INSTANCE.getLikelihood());
    }

    public double calcLogLikelihood() {
        return Math.log(SimpleNodeElimination.INSTANCE.getLikelihood());
    }

    public double calcFreeEnergy(int start, int end) {
        // only for debugging purposes
        if (MeanFieldForBayesNet.CALC_LIKELIHOOD) {
            return -calcLogLikelihood(start, end);
        }
        // #############################################

        // TODO: hier kann viel gecached werden, vor allem bei der Minimierung durch anpassen der Q
        double part1 = 0; // erste Summe der VE
        double part2 = 0; // zweite Summe der VE
        double tmp = 0;
        int hp; // speichert die durch l (condprob) kodierte Realisierung zum Elternindex p
        int hn = 0; // hiddenNodeNumber
        int hnp = 0; // hiddenNodeNumberParent

        BayesNetNode aktNode, aktParent;

        // Term 1, entropie
        for (int t = start; t <= end; t++) {
            for (int i = 0; i < bnh._myTreesArray[t].numberOfNodes; i++) {
                aktNode = bnh._myTreesArray[t].getNode(i);
                tmp = 0;
                if (!aktNode.props._isObserved) {
                    hn = node2hiddenNodeMap[aktNode.nodeNumber];
                    // laufe �ber alle Zust�nde in q[i]
                    for (int h = 0; h < 4; h++) {
                        if (q[hn + h] > 0) {
                            tmp += q[hn + h] * Math.log(q[hn + h]);
                        } else {
                            // System.err.println("1) log = -inf");
                        }
                    }
                }
                part1 += tmp;
            }
        }

        for (int t = start; t <= end; t++) {
            for (int i = 0; i < bnh._myTreesArray[t].numberOfNodes; i++) {
                aktNode = bnh._myTreesArray[t].getNode(i);
                int cpfSize = aktNode.CPF._size;
                double[][] PaktNodeCondProb;
                PaktNodeCondProb = aktNode.CPF.getCondProb(); // for fast access;
                // Util.visualizeMatrix(PaktNodeCondProb);
                // System.out.println("-------");
                int[][] PaktLookup = aktNode.CPF.getLookUpPointer();
                // System.out.println("aktNode: " + aktNode.getName());
                tmp = 0;
                hn = node2hiddenNodeMap[aktNode.nodeNumber];

                if (aktNode.isRoot() && !aktNode.props._isObserved) {
                    for (int h = 0; h < 4; h++) {
                        // System.out.println(q[hn][h] + " * log( " + aktNode.CPF.get(0, h) + " )");
                        if (PaktNodeCondProb[0][h] > 0) {// (aktNode.CPF.get(0, h) > 0)
                            tmp += q[hn + h] * Math.log(PaktNodeCondProb[0][h]); // (aktNode.CPF.get(0, h));
                        } else {
                            // tmp = Double.NEGATIVE_INFINITY;
                            // System.err.println("log = -inf");
                        }
                    }
                } else if (aktNode.props._isObserved) {
                    // laufe �ber condprob-array, l codiert die Beobachtung der eltern
                    // System.out.println("Fall 2 " + aktNode.getName());
                    for (int l = 0; l < cpfSize; l++) {
                        double prod = 1;
                        // int[] lookupParents = aktNode.CPF.getQueryByIndex(l);
                        for (int p = 0; p < aktNode.numberOfParents; p++) {
                            aktParent = aktNode.Aparents[p];
                            // System.out.println("parent: " + aktParent.getName());

                            hnp = node2hiddenNodeMap[aktParent.nodeNumber];
                            hp = PaktLookup[l][p]; // lookupParents[p];

                            if (!aktParent.props._isObserved) {
                                // das Produkt der q's berechnen
                                prod *= q[hnp + hp];
                            } else if (hp != aktParent.props._observedIndex) {
                                // prod *= aktParent.props._observation[hp]; //TODO: - is observable
                                prod = 0;
                            }
                        }
                        // System.out.println(l + ", " + aktNode.props._observedIndex+ ": " + prod + " * log( " +
                        // (aktNode.CPF.get(l,aktNode.props._observedIndex)) + " )");
                        if (PaktNodeCondProb[l][aktNode.props._observedIndex] > 0) { // (aktNode.CPF.get(l,
                                                                                     // aktNode.props._observedIndex) >
                                                                                     // 0) {
                            prod *= Math.log(PaktNodeCondProb[l][aktNode.props._observedIndex]);// aktNode.CPF.get(l,
                                                                                                // aktNode.props._observedIndex));
                        } else {
                            // System.err.println("2) log = -inf( "+ aktNode.CPF.get(l,aktNode.props._observedIndex)
                            // + " )");
                            prod = Double.NEGATIVE_INFINITY;
                        }
                        tmp += prod;
                    }
                } else /* if (!aktNode.props.isPhyloLeaf()) */{
                    for (int h2 = 0; h2 < 4; h2++) {
                        for (int l = 0; l < cpfSize; l++) {
                            double prod = 1;
                            // int[] lookupParents = aktNode.CPF.getQueryByIndex(l); // Merke Beobachtung zu l in den
                            // Eltern
                            for (int p = 0; p < aktNode.numberOfParents; p++) { // laufe �ber Eltern
                                aktParent = aktNode.Aparents[p]; // merke aktuellen Eltern knoten
                                hnp = node2hiddenNodeMap[aktParent.nodeNumber]; // merke hiddenNodeNumber
                                                                                // von den Eltern
                                hp = PaktLookup[l][p]; // lookupParents[p]; // Beobachtung des Eltern p merken
                                if (!aktParent.props._isObserved) {
                                    prod *= q[hnp + hp]; // das Produkt der q's berechnen: hnp = hiddennodenumber vom
                                                         // eltern, hp = Beobachtung in Eltern
                                } else if (hp != aktParent.props._observedIndex) {
                                    // prod *= aktParent.props._observation[hp]; //TODO: - is observable
                                    // TODO: Check this case
                                    prod = 0;
                                }
                            }
                            // System.out.println(prod + " * " + q[hn][h2] + " * log( " + (aktNode.CPF.get(l, h2)) +
                            // " )");
                            if (PaktNodeCondProb[l][h2] > 0) { // aktNode.CPF.get(l, h2)
                                prod *= q[hn + h2] * Math.log(PaktNodeCondProb[l][h2]); // aktNode.CPF.get(l, h2));
                            } else {
                                // System.err.println("3) log = -inf log( "+ aktNode.CPF.get(l, h2) + " )");
                                prod = Double.NEGATIVE_INFINITY;
                            }
                            tmp += prod;
                        }
                    }
                }
                // System.out.println(tmp);
                part2 += tmp;
            }
        }

        return part1 - part2;
    }
}
