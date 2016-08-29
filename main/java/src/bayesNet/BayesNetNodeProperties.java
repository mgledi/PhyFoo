package bayesNet;

import java.util.Arrays;
import java.util.logging.Logger;

import io.Alphabet;

public class BayesNetNodeProperties {
    public enum NODE_TYPE {
        PHYLO_LEAF_NODE, PHYLO_ROOT_NODE
    }

    private static Logger LOGGER = Logger.getLogger(BayesNetNodeProperties.class.getName());

    /**
     * Represents the type of the node.
     * 
     * @see NODE_TYPE
     */
    private NODE_TYPE nodeType;

    /** Parameter f�r MeanFieldForBayesNet-Caching */
    private int updateVersion = 1;

    /** vereinfacht �berpr�fung ob Knoten wirklich beobachtet wurde */
    public boolean _isObserved = false;

    /** repr�sentiert den index, der Beobachtung */
    public int _observedIndex;

    /** repr�sentiert vollst�ndige Beobachtung */
    public double[] _observation = new double[Alphabet.size];

    /** TRUE wenn Knoten ein Blatt eines Phylogenetischen Baums ist */
    public boolean isPhyloLeaf() {
        return nodeType == NODE_TYPE.PHYLO_LEAF_NODE;
    }

    public void setPhyloLeaf() {
        nodeType = NODE_TYPE.PHYLO_LEAF_NODE;
    }

    /** TRUE wenn Knoten eine Wurzel eines Phylogenetischen Baums ist */
    public boolean isPhyloRoot() {
        return nodeType == NODE_TYPE.PHYLO_ROOT_NODE;
    }

    public void setPhyloRoot() {
        nodeType = NODE_TYPE.PHYLO_ROOT_NODE;
    }

    /** TRUE wenn Knoten ein innerer Knoten eines Phylogenetischen Baums ist */
    public boolean isPhyloInnerNode() {
        return (nodeType != NODE_TYPE.PHYLO_LEAF_NODE && nodeType != NODE_TYPE.PHYLO_ROOT_NODE);
    }

    /** TRUE, wenn Knoten beobachtet wurde */
    public boolean isObserved() {
        return _isObserved;
    }

    /** gibt die Beobachtung des Knotens zur�ck */
    public double[] getObservation() {
        return _observation;
    }

    /** gibt den Alphabetindex der Beobachtung im Knoten zur�ck */
    public int getObservedIndex() {
        return _observedIndex;
    }

    /** setzt einen Knoten auf beobachtet */
    public void setObserved() {
        _isObserved = true;
    }

    /** setzt einen Knoten auf unbeobachtet */
    public void setUnobserved() {
        _isObserved = false;
    }

    /**
     * sets the given index
     * 
     * @param idx
     *            the index to set
     */
    public void setObservation(int idx) {
        // System.out.println(i + " " + Alphabet.size);
        Arrays.fill(_observation, 0);
        if (idx >= Alphabet.size) {
            // _myConsole.severe("Zusetzende Beobachtung au�erhalb des Alphabets");
            this.setUnobserved();
        } else {
            _observedIndex = idx;
            _observation[idx] = 1;
            this.setObserved();
        }
    }
}
