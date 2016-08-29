package bayesNet;

import java.util.ArrayList;

public class BayesNetNode {
    
    /** verwaltet die bedingten Wahrscheinlichkeiten */
    public CPF CPF;
    
    /**
     * enth�lt alle nicht generellen Eigenschaften f�r den BayesNetNode Bsp:
     * PhyloNode Observation Indexe, evtl. algorithmische Beziehungen
     */
    public BayesNetNodeProperties props = new BayesNetNodeProperties();

    /** Kinder des Knotens */
    public ArrayList<BayesNetNode> children = new ArrayList<BayesNetNode>(0);
    public BayesNetNode[] Achildren = new BayesNetNode[0]; // for fast access
    
    /** Eltern des Knotens */
    public ArrayList<BayesNetNode> parents = new ArrayList<BayesNetNode>(0);
    public BayesNetNode[] Aparents = new BayesNetNode[0]; // for fast access
    
    
    /** enth�lt die Anzahl der Kinder */
    public int numberOfChilds = 0;
    
    /** enth�lt die Anzahl der Eltern */
    public int numberOfParents = 0;
    
    /** nummer des knotens um direkt zuzugreifen */
    public int nodeNumber;
    
    /** Name des Knotens */
    protected String _name;
    
    /** Name des Knotens */
    public String _realname;
    
    /** Netz in dem der Knoten enthalten ist */
    protected BayesNet _myNet;

    /** enth�lt die Tiefe des Knotens */
    int depth;

    /** wird nur initialisiert, wenn eine der Knoten ein PhyloNode ist */
    protected double distanceToParent;
    

    /** Konstruktor ohne Namen */
    public BayesNetNode(BayesNet n) {
        init("noname" + n.nodes.size(), n);
    }

    /** Konstruktor mit Namen */
    public BayesNetNode(String name, BayesNet n) {
        init(name, n);
    }

    /** init Funktion, wird durch Konstruktor aufgerufen */
    private final void init(String name, BayesNet n) {
        // sicherheitshalber �berpr�fen ob der Knoten wirklich neu ist, das
        // hei�t die variablen noch nicht initialisiert wurden
        _name = name;
        _myNet = n;
        this.depth = 0;
        _myNet.addNode(this);
        this.CPF = new CPF(this);
    }

    /** holt die Tiefe des Knotens */
    public int getDepth() {
        return this.depth;
    }

    /** gibt den Name des Knotens zur�ck */
    public String getName() {
        return _name;
    }

    /** setzt den Namen des Knotens */
    public void setName(String name) {
        _name = name;
    }

    /** gibt den echten Namen des Knotens zur�ck */
    public String getRealName() {
        return _realname;
    }

    /** setzt den echten Namen des Knotens */
    public void setRealName(String name) {
        _realname = name;
    }
    
    /** setzt die entfernung zum elternknoten (nur bei Phylonode) */
    public void setDistanceToParent(double val) {
        // TODO: �berpr�fen ob eine Distanz gesetzt werden darf
        this.distanceToParent = val;
    }

    /** gibt die Distanz zum Elternknoten zur�ck */
    public double getDistanceToParent() {
        return this.distanceToParent;
    }

    /** f�gt dem Knoten ein Kind hinzu */
    public void addChild(BayesNetNode child) {
        this.children.add(child);
        this.numberOfChilds = this.children.size();
        child.parents.add(this);
        child.numberOfParents = child.parents.size();
        child.CPF.reinit();
        
        Achildren = children.toArray(new BayesNetNode[children.size()]);
        child.Aparents = child.parents.toArray(new BayesNetNode[child.parents.size()]);
    }

    /** entfernt ein kind von dem Knoten */
    public void removeChild(BayesNetNode child) {
        this.children.remove(child);
        this.numberOfChilds = this.children.size();
        child.parents.remove(this);
        child.numberOfParents = child.parents.size();
        child.CPF.reinit();
        
        Achildren = children.toArray(new BayesNetNode[children.size()]);
        child.Aparents = child.parents.toArray(new BayesNetNode[child.parents.size()]);
    }

    /** gibt ein Eltern nach seiner relativen Nummer zur�ck */
    public BayesNetNode getParent(int p) {
        return this.parents.get(p);
//        return Aparents[p];
    }

    /** gibt das Kind nach seiner relativen nummer zur�ck */
    public BayesNetNode getChild(int c) {
        return this.children.get(c);
//        return this.Achildren[c];
    }

    /** gibt true zur�ck, wenn aktueller Knoten eine Wurzel ist */
    public boolean isRoot() {
        return this.parents == null || this.parents.size() == 0;
    }

    /** gibt true zur�ck, wenn aktueller Knoten ein Kind ist */
    public boolean isLeaf() {
        return this.children == null || this.children.size() == 0;
    }

    /** gibt den index des �bergebenen Eltern zur�ck */
    public int getParentIndex(BayesNetNode parent) {
        return this.parents.indexOf(parent);
    }
    @Override
    public String toString() {
        return _name;
    }
}