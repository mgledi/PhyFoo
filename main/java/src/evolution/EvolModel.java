package evolution;

public abstract class EvolModel {
	/*
	 * The class of the Evolutionary model to use. TODO: don't use a static class
	 */
	public static String MODEL_CLASS = "FS81alpha";

	// TODO: must be configurable.
	protected int _alphabetSize = 4;


    /** @return the dimension of thie evolutionary model */
    public abstract int getDimension();

    /**
     * Reinitializes the transition matrix with the underlying stationary distributions and the given meta parameters,
     * e.g. the branch length
     * 
     * @param metaparam
     */
    public abstract void reinit(double... metaparam);

    /** @return a pointer to the stationary distributions */
    public abstract double[][] getStatDistr();

    /** Sets the given stationary distributions */
    public abstract void setStatDistr(double[]... A);

    /** Copies this EvolModel to an XML-encoded StringBuffer */
    public abstract StringBuffer toXML();

    public abstract double[][] getCurrentTransitions();

    @Deprecated
    public abstract double[][] getUnweightedTransitions();

    /**
     * This factory method creates a concrete instance of an {@link EvolModel}.
     * 
     * @param dimension
     * @return a concrete instance of an {@link EvolModel}
     */
    public static EvolModel getInstance(int dimension) {
		if (MODEL_CLASS.equals("FS81alpha")) {
            return new FS81alpha(dimension);
		} else if (MODEL_CLASS.equals("FS81beta")) {
            return new FS81beta(dimension);
		} else if (MODEL_CLASS.startsWith("HKY")) {
			double alpha_beta_ratio = Double.valueOf(MODEL_CLASS.split(":")[1]);
			return new HKY(dimension, alpha_beta_ratio);
        } else {
			throw new IllegalArgumentException("Can't handle class " + MODEL_CLASS);
        }
    }
}
