package bayesNet;

import java.util.ArrayList;
import java.util.logging.Logger;

import io.Alphabet;
import util.Util;

/**
 * CPF = Conditional Probability Function
 */
public class CPF {
	private Logger _logger = Logger.getLogger(this.getClass().getName());
	public static boolean ALLOW_RETURN_INDEPENDENT_TRANSITION = true;
	/** enth�lt einen Zeiger auf den besitzenden BayesKnoten */
	private BayesNetNode _myNode;

	/**
	 * enth�lt die Bedingten WKs, zeilen entsprechen den Beobachtungen in _myNode
	 */
	private double[][] _condProb;

	public double[][] _indepedentCondProb;

	/** enth�lt die Verbundwahrscheinlichkeit mit den Eltern */
	private double[][] _combinedProb;

	/**
	 * enth�lt die Gr��e der _condProb, abh�ngig vom Alphabet und der Menge der Eltern
	 */
	public int _size;

	/** enth�lt eine LookupTable, l kodiert Beobachtung der eltern */
	public int[][] _lookup;

	/**
	 * sizebuffer erleichtert die Konvertierung von der wahren position im Condprob-Array zu Beobachtung der eltern und umgedreht
	 */
	protected int[] _SizeBuffer;

	public CPF(BayesNetNode myNode) {
		this._myNode = myNode;
		reinit();
	}

	/** reinitialisiert die CPF, ACHTUNG: evtl enthaltene Werte gehen verloren */
	public void reinit() {
		_size = (int) Math.pow(Alphabet.size, _myNode.numberOfParents);
		_condProb = new double[_size][Alphabet.size];
		_combinedProb = new double[_size][Alphabet.size];
		_condProb = Util.getRandomStochMatrix(_size, Alphabet.size);

		// _SizeBuffer erzeugen ==> beschleunigt zugriff auf CondProb array
		_SizeBuffer = new int[_myNode.numberOfParents];
		int size = 1;
		for (int i = 0; i < _myNode.numberOfParents; i++) {
			_SizeBuffer[i] = size;
			size *= Alphabet.size;
		}

		// Lookuptable erzeugen, um sp�ter auf Realisierungen der Eltern per
		// Index zugreifen zu k�nnen
		this.generateLookUpTable();
	}

	/**
	 * berechnet die Verbundwahrscheinlichkeit des aktuellen Knotens mit seinen Eltern
	 */
	public void calculateCombinedProb() {
		double[] parentMargin;
		int[] query;
		BayesNetNode parent;
		_combinedProb = Util.arraycopy(_condProb);
		// System.out.println("---------");
		// Util.visualizeMatrix(_combinedProb, 4, true);
		if (_myNode.numberOfParents == 0) {
			// nothing
		} else {
			// Eltern m�ssen aufmultipliziert werden
			for (int l = 0; l < _size; l++) {
				query = getQueryByIndex(l);
				for (int p = 0; p < _myNode.numberOfParents; p++) {
					parent = _myNode.getParent(p);
					parentMargin = parent.CPF.getFullMarginalizedDistr();
					for (int a = 0; a < Alphabet.size; a++) {
						// System.out.println("P(X="+a+"|Y="+query[p]+") * P(Y="+query[p]+")");
						_combinedProb[l][a] *= parentMargin[query[p]];
					}
				}
			}
		}
		// Util.visualizeMatrix(_combinedProb, 12, true);
	}

	/** */
	public double[] getFullMarginalizedDistr() {
		double[] margin = new double[Alphabet.size];
		for (int l1 = 0; l1 < _size; l1++) {
			for (int a = 0; a < Alphabet.size; a++) {
				margin[a] += _combinedProb[l1][a];
			}
		}
		// System.out.println(Util.array2string(" ,", margin));
		return margin;
	}

	/** berechnet die Randverteilung �ber die Eltern */
	public void calculateMarginalDistr() {
		this.calculateMarginDistr(_myNode.parents);
	}

	/** berechnet die Randverteilung der �bergebenen Knoten */
	private void calculateMarginDistr(ArrayList<BayesNetNode> nodes) {
		if (nodes.indexOf(_myNode) != -1) {
			_logger.severe("Bei der Berechnung der Randverteilung ist ein Fehler aufgetreten!");
		}

		int size = (int) (_size / Math.pow(Alphabet.size, nodes.size()));
		double[][] margin = new double[size][Alphabet.size];

		// laufe �ber die ganze Verbundwahrscheinlichkeit und summiere
		// entsprechend auf
		for (int l1 = 0; l1 < _size; l1++) {
			for (int a = 0; a < Alphabet.size; a++) {
				margin[0][a] += _combinedProb[l1][a];
				// TODO:
				// for(int i=0; i < nodes.size(); i++) {
				// pos = _myNode.getParentIndex(nodes.get(i));
				// if(pos >= 0) {
				//
				//
				// }
				// }
			}
		}
		// System.out.println(Arrays.toString(margin[0]));

	}

	/** gibt die VerbundWK zur�ck, diese sollte vorher berechnet werden */
	public double[][] getCombinedProb() {
		return _combinedProb;
	}

	/** returns a pointer to the conditional probabilty matrix */
	public double[][] getCondProb() {
		if (_myNode.props.isPhyloLeaf() && !_myNode.props._isObserved && ALLOW_RETURN_INDEPENDENT_TRANSITION) {
			return _indepedentCondProb;
		} else {
			return _condProb;
		}
	}

	/** setzt neue bedingte Wahrscheinlichkeiten f�r die CPF */
	public void setCondProb(double[][] condProb) {
		// �berpr�fen ob die zu setzende bedingte WK, die gleiche Dimension hat
		if (condProb.length != _condProb.length || condProb[0].length != _condProb[0].length) {
			System.out.println(condProb.length + " : " + condProb[0].length + " expected " + _condProb.length + " : " + _condProb[0].length);
			_logger.severe("Falsche Dimension in setCondProb f�r Knoten " + _myNode.getName());
		} else {
			_condProb = condProb;
		}
	}

	/** gibt die vorberechnete Gr��e der CPF zur�ck, (Alphabet ^ Eltern) */
	public int size() {
		return _size;
	}

	/**
	 * gibt den _sizeBuffer zur�ck, welcher einen schnelleren Sprung durch das CondProb array erm�glicht
	 */
	public int[] getSizeBuffer() {
		return _SizeBuffer;
	}

	protected int getIndexByQuery(int[] query) {
		if (query.length != _myNode.numberOfParents) {
			_logger.severe("get() Falsche Dimension der Anfrage = (" + query.length + "). " + "erwartet: " + _myNode.numberOfParents);
			return -1;
		} else {
			int l = 0;
			for (int i = 0; i < query.length; i++) {
				l += _SizeBuffer[i] * query[i];
			}
			return l;
		}
	}

	/**
	 * gibt den zugeh�rigen Wert aus dem CondProb array zur�ck, query verschl�sselt die Adresse
	 */
	public double get(int[] query, int obs) {
		int l = getIndexByQuery(query);
		if (l == -1) {
			_logger.severe("get(): Es konnte keine Position im Array bestimmt werden");
			return -1;
		} else {
			return _condProb[l][obs];
		}
	}

	/** gibt den zugeh�rigen Vektor aus dem CondProb array zur�ck */
	public double[] get(int[] query) {
		int l = getIndexByQuery(query);
		if (l == -1) {
			_logger.severe("get(): Es konnte keine Position im Array bestimmt werden");
			return null;
		} else {

		}
		return _condProb[l];
	}

	/**
	 * gibt den zugeh�rigen Wert aus dem CondProb array zur�ck, l ist ihr die direkte Adresse
	 */
	public double get(int l, int obs) {
		if (_myNode.props.isPhyloLeaf() && !_myNode.props._isObserved && ALLOW_RETURN_INDEPENDENT_TRANSITION) {
			return _indepedentCondProb[l][obs];
		} else {
			return _condProb[l][obs];
		}
	}

	public int[][] getLookUpPointer() {
		return _lookup;
	}

	public int[] getQueryByIndex(int l) {
		return this.getQueryByIndex(l, true);
	}

	/**
	 * erzeugt aus einem Index ein Beobachtungsarray der Eltern dies kann w�hrend der Laufzeit geschehen oder bereits vor alle m�glichen Indices vorberechnet
	 * werden. Aufgrund der Modulo-Operation sehr teuer.
	 */
	public int[] getQueryByIndex(int l, boolean cache) {
		// wenn cache = true, dann gib den vorkalkulierten wert zur�ck
		if (cache) {
			return _lookup[l];
		}

		// initialisiere Anfrage
		int[] query = new int[_myNode.numberOfParents];

		// laufe �ber all Eltern um aus dem �bergebenen Index die beobachtung zu
		// berechnung
		for (int i = _myNode.numberOfParents - 1; i >= 0; i--) {
			query[i] = (int) l / _SizeBuffer[i] % Alphabet.size; // TODO: testen
			if (i == 0) {
				query[0] = (int) l % Alphabet.size;
			}
		}
		return query;
	}

	/** �berpr�ft ob die �bergebene beobachtung der Wahrheit entspricht */
	public boolean isTrueObs(int l) {
		int[] query = getQueryByIndex(l);
		for (int p = 0; p < _myNode.numberOfParents; p++) {
			if (_myNode.getParent(p).props.isObserved() && _myNode.getParent(p).props.getObservedIndex() == query[p]) {
				return false;
			}
		}
		return true;
	}

	/** �berpr�ft ob die �bergebene beobachtung der Wahrheit entspricht */
	public boolean isTrueObs(int[] query) {
		for (int p = 0; p < _myNode.numberOfParents; p++) {
			if (_myNode.getParent(p).props.isObserved() && _myNode.getParent(p).props.getObservedIndex() == query[p]) {
				return false;
			}
		}
		return true;
	}

	/** erzeugt eine Lookuptable, l kodiert dabei eine Beobachtung in den eltern */
	private void generateLookUpTable() {
		// initialisiere LookupTable
		_lookup = new int[_size][_myNode.numberOfParents];
		// laufe �ber alle m�glichen indexe und generiere zu jedem die
		// Beobachtung
		for (int l = 0; l < _size; l++) {
			_lookup[l] = this.getQueryByIndex(l, false);
		}
		// Util.visualizeMatrix(_lookup, 1, true);
	}
}
