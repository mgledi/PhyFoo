package bayesNet;

import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bayesNet.BayesNet;
import bayesNet.BayesNetNode;
import bayesNet.CPF;
import optimizing.TemperatureBasedMotifParamOptimizer;
import util.Util;
import de.jstacs.utils.DoubleList;
import de.jstacs.utils.Normalisation;
import evolution.EvolModel;

/** Virtual environment in a BayesNet. Represents all nodes within a phylogenetic tree */
public class VirtualTree {

    /**
     * This variable indicates a filtering process on species 1. It is only used in
     * {@link TemperatureBasedMotifParamOptimizer}.
     */
    public double TEMPERATURE = 1;

    /** A Vector of all nodes in this tree */
    private ArrayList<BayesNetNode> myNodes = new ArrayList<BayesNetNode>(0);

    /** A vector with pointers to all leafs within this tree */
    private ArrayList<BayesNetNode> myLeafs = new ArrayList<BayesNetNode>(0);

    private String newickString;

    /** the underlying evolutionary model of this tree. */
    private EvolModel evolutionaryModel;

    /** number of nodes in this tree */
    public int numberOfNodes;

    /** number of leafs in this tree */
    public int numberOfLeafs;

    /** A pointer to the parent container. */
    BayesNet net;

    public VirtualTree(BayesNet net, String newickString) {
        this.net = net;
        for (int i = 0; i < net.numberOfNodes; i++) {
            this.myNodes.add(net.getNode(i)); // ACHTUNG: nur Referenzen sollen gespeichert werden
            if (net.getNode(i).isLeaf()) {
                this.myLeafs.add(net.getNode(i));
            }
        }
        this.numberOfNodes = this.myNodes.size();
        this.numberOfLeafs = this.myLeafs.size();
        this.newickString = newickString;
    }

    /**
     * If not set, determines the newick string from the tree and returns it.
     * 
     * @return the newick string
     */
    public String getNewickString() {
        newickString = this.generateNewickString();
        return newickString;
    }

    /** returns a pointer to the used evolutionary model */
    public EvolModel getEvolModel() {
        return evolutionaryModel;
    }

    /** sets the evolutionary model */
    public void setEvolModel(EvolModel gfs) {
        evolutionaryModel = gfs;
    }

    /**
     * Copies all relevant parameters from the underlying evolutionary model to the conditional probability function
     * {@link CPF}
     */
    public void initParametersFromGF() {
        BayesNetNode aktNode;
        // if TEMPERATURE is 1, then no filtering is assumed
        if (TEMPERATURE == 1) {
            for (int i = 0; i < numberOfNodes; i++) {
                aktNode = myNodes.get(i);
                if (aktNode.props.isPhyloRoot()) {
                    aktNode.CPF.setCondProb(Util.arraycopy(evolutionaryModel.getStatDistr()));
                    aktNode.CPF._indepedentCondProb = Util.arraycopy(evolutionaryModel.getStatDistr());
                } else {
                    evolutionaryModel.reinit(aktNode.getDistanceToParent());
                    aktNode.CPF.setCondProb(Util.arraycopy(evolutionaryModel.getCurrentTransitions()));
                    aktNode.CPF._indepedentCondProb = evolutionaryModel.getUnweightedTransitions();
                }
            }
        } else {
            // TODO fix this
            if (evolutionaryModel.getStatDistr().length > 1) {
                throw new IllegalArgumentException(
                        "Temperature based parameter setting only allowed for Markov Models of order 0.");
            }
            if (numberOfLeafs +1 > numberOfNodes) {
                throw new IllegalArgumentException(
                        "A star topology is expected. In a star topology there exists only one non leaf node, the root");
            }

            // ###################### Step 1)
            // initialize conditional probabilities for species 1
            double[][] pi = evolutionaryModel.getStatDistr();
            double[][] pi_A = new double[pi.length][];
            // calculate the new sharpened vector
            for (int k = 0; k < pi_A.length; k++) {
                pi_A[k] = new double[pi[k].length];
                for (int a = 0; a < pi_A[k].length; a++) {
                    pi_A[k][a] = Math.pow(pi[k][a], TEMPERATURE);
                }
                Normalisation.sumNormalisation(pi_A[k]);
            }

            {   // TODO: for all nodes up to the root
                // TODO: respect order > 0
                // ###################### Step 2)
                // mutate distribution in species 1 to obtain distribution for the root (premordial species)
                evolutionaryModel.reinit(this.getLeaf(0).getDistanceToParent()); // use unfiltered mutation matrix
                double[][] F_pi = evolutionaryModel.getCurrentTransitions();
                double[][] pi_Y = new double[1][];
                pi_Y[0] = times(pi_A[0], F_pi);
                this.getNode(0).CPF.setCondProb(pi_Y);
    
                // ###################### Step 3)
                // determine P(A|Y) from pi_A, pi_Y, and F_pi
                double[][] reverse_F_pi = Util.getReverseFS81(pi_A[0], pi_Y[0], F_pi);
                this.getLeaf(0).CPF.setCondProb(reverse_F_pi);
                
                for(int o = 1; o < this.numberOfLeafs; o++) {
                    evolutionaryModel.reinit(getLeaf(o).getDistanceToParent());
                    getLeaf(o).CPF.setCondProb(Util.arraycopy(evolutionaryModel.getCurrentTransitions()));
                }
            }
        }
    }
    
    public static double[] times(double[] vector, double[][] matrix) {
        double[] result = new double[vector.length];
        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < matrix[0].length; col++) {
                result[col] += vector[row] * matrix[row][col];
            }
        }
        return result;
    }

    /** @return node with index i */
    public BayesNetNode getNode(int i) {
        return myNodes.get(i);
    }

    /** @return leaf with index i */
    public BayesNetNode getLeaf(int i) {
        return myLeafs.get(i);
    }

    /** @return true if the given {@link BayesNetNode} exists */
    public boolean contains(BayesNetNode bn) {
        return myNodes.contains(bn);
    }

    /** @return index of the given {@link BayesNetNode} */
    public int getNodePos(BayesNetNode bn) {
        return myNodes.indexOf(bn);
    }

    /** gibt die Knoten als Liste aus */
    public void outNodeList() {
        for (int k = 0; k < this.numberOfNodes; k++) {
            System.out.println(this.getNode(k).getName() + " (" + this.getNode(k).nodeNumber + ") " + " Parents: "
                    + this.getNode(k).numberOfParents + " Childs: " + this.getNode(k).numberOfChilds + " ");
        }
    }

    /** setzt die �bergebene Beobachtung in die Bl�tter */
    public void setObservation(int[] obs) {
        // this should be the standard case
        if(obs.length == myLeafs.size()) {
//            System.out.println("Set Observation for leaves: ");
            for (int i = 0; i < myLeafs.size(); i++) {
//                System.out.println("i " + obs[i]);
                myLeafs.get(i).props.setObservation(obs[i]);
            }
        } else if(obs.length == myNodes.size()) {
            for (int i = 0; i < myNodes.size(); i++) {
                myNodes.get(i).props.setObservation(obs[i]);
            }
        } else {
            throw new IllegalArgumentException("Can not map Observation to tree.");
        }
    }

    public int[] getObservation() {
        int[] obs = new int[this.numberOfLeafs];
        for (int i = 0; i < this.numberOfLeafs; i++) {
            obs[i] = myLeafs.get(i).props.getObservedIndex();
        }
        return obs;
    }

    public String generateNewickString() {
        return this.generateNewickString(this.getNode(0), 0) + ";";
    }

    public String generateNewickString(BayesNetNode aktnode, int depth) {
        boolean hasChildren = false;
        int rc = 0;
        StringBuilder sb = new StringBuilder();
        for (int c = 0; c < aktnode.numberOfChilds; c++) {
            BayesNetNode child = aktnode.getChild(c);
            if (myNodes.contains(child)) {
                hasChildren = true;
                if (rc > 0) {
                    sb.append(",");
                }
                sb.append(this.generateNewickString(aktnode.getChild(c), depth + 1)).append(":")
                        .append(Util.round(child.getDistanceToParent(), 3));
                rc++;
            }
        }
        if (!hasChildren) {
            return aktnode.getRealName();
        }
        sb.insert(0, "(");
        sb.append(")");
        return sb.toString();
    }

    /** Reinits for the egdelenghts of the underlying phylogenetic trees. */
    public void reinitEdglengths(String newick) {
        // TODO: check for each tree correctness of topology against actual topology
        DoubleList dl = new DoubleList();
        Pattern p = Pattern.compile(":([0-9]*\\.[0-9]*)");
        Matcher m = p.matcher(newick);
        int i = 0;
        while (m.find()) {
            dl.add(Double.valueOf(m.group(1)));
        }

        i = 0;
        for (int n = 0; n < numberOfNodes; n++) {
            BayesNetNode node = getNode(n);
            if (!node.props.isPhyloRoot()) {
                node.setDistanceToParent(dl.get(i++));
            }
        }
        this.initParametersFromGF();
    }
}
