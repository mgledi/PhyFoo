package bayesNet;

import java.util.ArrayList;
import java.util.logging.Logger;

/**
 * 
 * @author m.nettling
 * 
 */
public class BayesNet {
	private static Logger logger = Logger.getLogger(BayesNet.class.getName());

	/** The name of the bayes net */
	private String name;

	/** A pointer to the root node of the bayesnet */
	private BayesNetNode root;

	/** An array of all nodes in this BayesNet */
	public ArrayList<BayesNetNode> nodes = new ArrayList<BayesNetNode>(0);

	/** The number of nodes */
	public int numberOfNodes = 0;

	public BayesNet() {
		this.init(null, "BayesNet");
	}

	public BayesNet(String name) {
		this.init(null, name);
	}

	public BayesNet(BayesNetNode root, String name) {
		this.init(root, name);
	}

	public BayesNet(BayesNetNode root) {
		this.init(root, "BayesNet");
	}

	protected void init(BayesNetNode root, String name) {
		if (root != null) {
			this.setRoot(root);
		}
		this.name = name;
	}

	public void setRoot(BayesNetNode root) {
		if (this.root != null) {
		} else {
			this.addNode(root);
			this.root = root;
		}
	}

	public BayesNetNode getRoot() {
		return this.root;
	}

	public String getName() {
		return this.name;
	}

	public int getNumberOfNodes() {
		return this.numberOfNodes;
	}

	public BayesNetNode getNode(int i) {
		return this.nodes.get(i);
	}

	public boolean checkNodeExists(BayesNetNode node) {
		if (nodes.contains(node)) {
			return true;
		}

		for (int i = 0; i < this.numberOfNodes; i++) {
			if (node.getName().equals(this.nodes.get(i).getName())) {
				logger.warning("Ein anderes Objekt der Knotenliste (" + node.getName() + ") hat bereits diesen Namen.");
				return true;
			}
		}
		return false;
	}

	public boolean addNode(BayesNetNode node) {
		if (!checkNodeExists(node)) {
			this.nodes.add(node);
			this.numberOfNodes = this.nodes.size();
			node.nodeNumber = this.numberOfNodes - 1;
			return true;
		} else {
			// logger.severe("Ein Knoten mit dem Namen " + node.getName() +
			// " exisitiert bereits und wird nicht hinzugef�gt.");
		}
		return false;
	}

	public BayesNetNode getNodeByName(String name) {
		for (int i = 0; i < this.nodes.size(); i++) {
			if (name.equals(this.nodes.get(i).getName())) {
				return this.nodes.get(i);
			}
		}
		return null;
	}

	public void outNodeList() {
		for (int k = 0; k < this.numberOfNodes; k++) {
			System.out.println(this.getNode(k).getName() + " (" + this.getNode(k).nodeNumber + ") " + " Parents: " + this.getNode(k).numberOfParents
			        + " Childs: " + this.getNode(k).numberOfChilds + " ");
		}
	}
}
