package bayesNet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import de.jstacs.NonParsableException;
import de.jstacs.io.XMLParser;
import evolution.FS81alpha;
import util.Util;

/**
 * Container for {@link BayesNet}, {@link BayesNetNode}s, and {@link VirtualTree}s.
 * 
 * @author Martin NEttling
 */
public class BayesNetHandler {
    private static Logger logger = Logger.getLogger(BayesNetHandler.class.getName());

    private BayesNet _myNet;
	private List<VirtualTree> _myTrees;
    public VirtualTree[] _myTreesArray;
    private int[][] _connectionTable;
    public int motifLength = 0;

	// ######### variables for topological sorting #####################
	private int[] maxNodeDepth;
    private int maxDepth = 0;

    public BayesNetHandler(BayesNet myNet) {
        _myNet = myNet;
		_myTrees = new ArrayList<>();
    }

    /**
     * Constructor: loads the Bayesnet and the BayesNetHandler from the given
     * XML
     * 
     * @throws NonParsableException
     */
    public BayesNetHandler(StringBuffer xml) throws NonParsableException {
        fromXML(xml);
    }

    /** f�gt dem aktuellen Bayesnetz ein SubNetz hinzu */
    public void addBayesNet(BayesNet net, String newickString) {
        for (int i = 0; i < net.numberOfNodes; i++) {
            _myNet.addNode(net.getNode(i));
        }
        motifLength++;
        _myTrees.add(new VirtualTree(net, newickString));
        _myTreesArray = _myTrees.toArray(new VirtualTree[0]);
        _connectionTable = this.generateNewConnectionTable();
    }

    /**
     * gibt eine Tabelle zur�ck, die die Verkn�pfung der Virtuellen b�ume
     * darstellt
     */
    public int[][] getConnectionTable() {
        return _connectionTable;
    }

    /** gibt alle Subgraphen (VirtualTrees) als Array zur�ck */
    public VirtualTree[] getVirtualTrees() {
        return _myTreesArray;
    }

    /** gibt ein Subnetz (virtualtree) zur�ck */
    public VirtualTree getVirtualTree(int i) {
        return _myTreesArray[i];
    }

    /** setzt den ersten knoten in der Knotenliste als Wurzel */
    public void setSimpleRoot() {
        for (int i = 0; i < _myNet.numberOfNodes; i++) {
            if (_myNet.getNode(i).numberOfParents == 0) {
                _myNet.setRoot(_myNet.getNode(0));
                return;
            }
        }
    }

    /** gibt einen Zeiger auf das Zugrunde liegende BayesNetz zur�ck */
    public BayesNet getNet() {
        return _myNet;
    }

    /** setzt an der motivposition pos die �bergebene Beobachtung obs */
    public void setObservation(int[] obs, int pos) {
        _myTrees.get(pos).setObservation(obs);
    }

    /** l�uft systematisch �ber das Netz und l�sst alle Verbundwks berechnen */
    public void calculateCombinedProbs() {
        // es wird von einer Topologischen Sortierung ausgegangen
        for (int i = 0; i < _myNet.numberOfNodes; i++) {
            _myNet.getNode(i).CPF.calculateCombinedProb();
        }
    }

    /** erzeugt zu den Paramteren eines Bayesnetzes eine volls�ndige Beobachtung */
    public int[] drawFullObservation() {
        int[] oldObs = new int[_myNet.numberOfNodes];
        Arrays.fill(oldObs, -1);
        return drawFullObservation(oldObs);
    }

    /**
     * generates a full observation for this net using internal probabilities. If an observation is set for a node at
     * index node_number in oldObs, this observation is used.
     */
    public int[] drawFullObservation(int[] oldObs) {
        int[] fullObs = new int[_myNet.numberOfNodes];
        Arrays.fill(fullObs, -1);

        // es wird von einer Topologischen Sortierung ausgegangen
        for (int i = 0; i < _myNet.numberOfNodes; i++) {
            // wenn bereits eine Beobachtung bekannt ist, dann diese verwenden
            if (oldObs[_myNet.getNode(i).nodeNumber] != -1) {
                fullObs[_myNet.getNode(i).nodeNumber] = oldObs[_myNet.getNode(i).nodeNumber];
            } else {
                fullObs[_myNet.getNode(i).nodeNumber] = drawObservationForNode(_myNet.getNode(i), fullObs);
            }
        }
        // System.out.println(Util.array2string(",", fullObs));

        return fullObs;
    }

    private int drawObservationForNode(BayesNetNode aktNode, int[] fullObs) {
        int[] query = new int[aktNode.numberOfParents];
        for (int p = 0; p < aktNode.numberOfParents; p++) {
            query[p] = fullObs[aktNode.getParent(p).nodeNumber];
        }

        double r = Util.getRandomDouble(0, 1);
        // System.out.println("Distr: ("+aktNode.nodeNumber+")" +
        // Arrays.toString(aktNode.CPF.get(query)));
        return Util.getRandomIndex(aktNode.CPF.get(query), r);
    }

    // ################################################################################################################
    // ####################################### Methoden zur topologischen
    // Sortierung der Knoten im Bayesnetz ########
    public void topologicalSort() {
        ArrayList<BayesNetNode> newNodeList = new ArrayList<BayesNetNode>(_myNet.numberOfNodes);
        this.calcMaxDepthArray();
        for (int d = 0; d <= this.maxDepth; d++) {
            for (int i = 0; i < _myNet.numberOfNodes; i++) {
                if (d == this.maxNodeDepth[_myNet.getNode(i).nodeNumber]) {
                    newNodeList.add(_myNet.getNode(i));
                }
            }
        }
        // synchronisiert Indexe und NodeNumbers, erlaubt schnelleren und
        // saubereren Zugriff
        for (int i = 0; i < _myNet.numberOfNodes; i++) {
            newNodeList.get(i).nodeNumber = i;
        }
        _myNet.nodes = newNodeList;
    }

    /**
     * berechnet f�r jeden Knoten die maximale Tiefe im Netz von allen Wurzeln
     * aus gesehen
     */
    private void calcMaxDepthArray() {
        this.maxNodeDepth = new int[_myNet.numberOfNodes];
        Arrays.fill(this.maxNodeDepth, 0);
        for (int i = 0; i < _myNet.numberOfNodes; i++) {
            if (_myNet.getNode(i).isRoot()) {
                this.calcMaxDepthArray(_myNet.getNode(i), 0);
            }
        }
    }

    /**
     * berechnet f�r jeden Knoten unter aktNode die Tiefe im Netz, ruft sich
     * rekursiv auf
     */
    private void calcMaxDepthArray(BayesNetNode aktNode, int depth) {
        if (depth > this.maxDepth) {
            this.maxDepth = depth;
        }

        // da ein Bayesnetz ein gerichteter Graph ist, muss keine �berpr�fung
        // auf eine Rekursionsschleife vorgenommen werden
        if (this.maxNodeDepth[aktNode.nodeNumber] <= depth) {
            this.maxNodeDepth[aktNode.nodeNumber] = depth;
            for (int i = 0; i < aktNode.numberOfChilds; i++) {
                this.calcMaxDepthArray(aktNode.getChild(i), depth + 1);
            }
        }
    }

    public void seperateVirtualTrees(int t1, int t2) {
        this.seperateVirtualTrees(_myTrees.get(t1), _myTrees.get(t2));
    }

    public void seperateVirtualTrees(VirtualTree t1, VirtualTree t2) {
        int pos1 = _myTrees.indexOf(t1);
        int pos2 = _myTrees.indexOf(t2);

        if (t1.numberOfNodes != t2.numberOfNodes) {
            logger.fine("seperateVirtualTrees: B�ume scheinen unterschiedliche Struktur zu haben");
            return;
        }

        for (int i = 0; i < t1.numberOfNodes; i++) {
            // horizontale Kante setzen
            t1.getNode(i).removeChild(t2.getNode(i));
        }
        _connectionTable[pos1][pos2] = 0;
    }

    /**
     * verkn�pft alle Konten der �bergebenen Virtual Trees paarweise horizontal
     * das hei�t: Knoten l von t2 wird Kind von Knoten l von t1
     */
    public void connectVirtualTrees(int t1, int t2) {
        this.connectVirtualTrees(_myTrees.get(t1), _myTrees.get(t2));
    }

    /**
     * verkn�pft alle Konten der �bergebenen Virtual Trees paarweise horizontal
     * das hei�t: Knoten l von t2 wird Kind von Knoten l von t1
     */
    public void connectVirtualTrees(VirtualTree t1, VirtualTree t2) {
        if (t1.numberOfNodes != t2.numberOfNodes) {
            logger.fine("connectVirtualTrees: B�ume scheinen unterschiedliche Struktur zu haben");
            return;
        }

        // merken welche B�ume bereits konkatiniert wurden
        int[][] cTable = generateNewConnectionTable();
        cTable[_myTrees.indexOf(t1)][_myTrees.indexOf(t2)] = 1;
        _connectionTable = cTable;

        if (checkForDimerConsistency(cTable)) {
        } else {
            logger.fine("connectVirtualTrees(): Dimerkonsitenz verletzt.");
        }

        for (int i = 0; i < t1.numberOfNodes; i++) {
            // horizontale Kante setzen
            t1.getNode(i).addChild(t2.getNode(i));

            // Diagonale Kante setzen
            if (!t2.getNode(i).isLeaf())
                for (int c = 0; c < t2.getNode(i).numberOfChilds; c++) {
                    // t1.getNode(i).addChild(t2.getNode(i).getChild(c));
                }
        }
    }

    /** �berpr�ft, ob durch das Modell wirklich nur Dimere behandelt werden */
    private boolean checkForDimerConsistency(int[][] cTable) {
        int[] ingoingEdges = new int[motifLength];
        for (int i = 0; i < motifLength; i++) {
            for (int j = 0; j < motifLength; j++) {
                ingoingEdges[j] += cTable[i][j];
                if (ingoingEdges[j] > 1) {
                    return false;
                }
            }
        }
        return true;
    }

    /**
     * generiert eine neue ConnectionTable, so dass die alte darin erhalten
     * bleibt
     */
    private int[][] generateNewConnectionTable() {
        int[][] newConnectionTable = new int[motifLength][motifLength];
        if (_connectionTable != null) {
            for (int i = 0; i < _connectionTable.length; i++) {
                for (int j = 0; j < _connectionTable.length; j++) {
                    newConnectionTable[i][j] = _connectionTable[i][j];
                }
            }
        }
        return newConnectionTable;
    }

    public StringBuffer toXML() {
        StringBuffer xml = new StringBuffer();
        StringBuffer Stmp;

        // Anzahl Positionen also VirtualTrees
        XMLParser.appendObjectWithTags(xml, motifLength, "numberOfTrees");

        // Abh�ngigkeit der B�ume speichern
        XMLParser.appendObjectWithTags(xml, _connectionTable, "connectionTable");

        // Virtuelle B�ume mit Eigenschaften speichern
        for (int i = 0; i < _myTrees.size(); i++) {
            Stmp = getVirtualTree(i).getEvolModel().toXML();
            XMLParser.addTags(Stmp, "evolModel");
            XMLParser.appendObjectWithTags(Stmp, getVirtualTree(i).TEMPERATURE, "TEMPERATURE");
            XMLParser.appendObjectWithTags(Stmp, getVirtualTree(i).getNewickString(), "newickString");
            XMLParser.addTags(Stmp, "VirtualTree" + i);
            xml.append(Stmp);
        }
        return xml;
    }

    /**
     * l�d aus einem XML-string alte parameter in das netz. ACHTUNG: die
     * Struktur muss mit der im XML-File �bereinstimmen
     * 
     * @param xml
     * @throws NonParsableException
     */
    public void fromXML(StringBuffer xml) throws NonParsableException {
        int ml = (Integer) XMLParser.extractObjectForTags(xml, "numberOfTrees");
        int[][] ct = (int[][]) XMLParser.extractObjectForTags(xml, "connectionTable");

        _myTrees = new ArrayList<VirtualTree>(ml);
        _myNet = new BayesNet(); // TODO: save also name
        // load temporary connectionTable
        StringBuffer Stmp;
        FS81alpha[] gfs = new FS81alpha[ml];
        double[] temperatures = new double[ml];
        // Load Basic Structure
        for (int i = 0; i < ml; i++) {
            Stmp = XMLParser.extractForTag(xml, "VirtualTree" + i);
            String newick = (String) XMLParser.extractObjectForTags(Stmp, "newickString");
            if (XMLParser.hasTag(Stmp, "TEMPERATURE", null, null)) {
                temperatures[i] = (Double) XMLParser.extractObjectForTags(Stmp, "TEMPERATURE");
            } else {
                temperatures[i] = 1;
            }
            this.addBayesNet(io.NewickToBayesNet.getTree(newick, "pos_" + i), newick);
            gfs[i] = FS81alpha.getObjectFromXML(XMLParser.extractForTag(Stmp, "evolModel"));
        }
        // System.out.println("trees loaded. Connect trees.");

        // connect trees vom stored connectionTable
        for (int i = ct.length - 1; i >= 0; i--) {
            for (int j = 0; j < ct[i].length; j++) {
                if (ct[i][j] == 1) {
                    this.connectVirtualTrees(i, j);
                }
            }
        }
        // System.out.println("trees connected");

        for (int i = 0; i < ml; i++) {
            VirtualTree vt = this.getVirtualTree(i);
            vt.TEMPERATURE = temperatures[i];
            vt.setEvolModel(gfs[i]);
            vt.initParametersFromGF();
        }
    }
}
