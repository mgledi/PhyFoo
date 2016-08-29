package util;

import java.io.ObjectInputStream.GetField;
import java.util.Properties;

public class Parameters {
    private Properties props;

    public enum ModelType {
        FOREGROUND("_FG"), BACKGROUND("_BG");

        protected final String parameterSuffix;

        ModelType(String parameterSuffix) {
            this.parameterSuffix = parameterSuffix;
        }
    }

    public Parameters(String[] args) {
        this(args, new Properties(), false);
    }

    /**
     * @param args
     *            parameters to overwrite in the given props object
     * @param props
     *            the properties-object to fill
     * @param argsParametersMustExist
     *            if true, the given properties must contain all keys also used in args
     */
    public Parameters(String[] args, Properties props, boolean argsParametersMustExist) {
        this.props = props;
        this.addParameters(args);
    }

    public void addParameters(String ... args) {
        Config.parseProperties(props, args, false);
    }

    public String getProperty(String key, String defaultValue) {
        return props.getProperty(key, defaultValue);
    }
    
    public int getPropertyAsInt(String key, int defaultValue) {
        return Integer.valueOf(props.getProperty(key, String.valueOf(defaultValue)));
    }
    
    public double getPropertyAsDouble(String key, double defaultValue) {
        return Double.valueOf(props.getProperty(key, String.valueOf(defaultValue)));
    }
    
    public String getProperty(String key, boolean allowNull) {
        String val = props.getProperty(key);
        if(val == null && allowNull) {
            return val;
        } else {
            throw new IllegalArgumentException("No value found for key " + key);
        }
    }
    
    public int getPropertyAsInt(String key, boolean allowNull) {
        String val = props.getProperty(key);
        if(val == null && allowNull) {
            return Integer.valueOf(val);
        } else {
            throw new IllegalArgumentException("No value found for key " + key);
        }
    }

    public double getPropertyAsDouble(String key, boolean allowNull) {
        String val = props.getProperty(key);
        if(val == null && allowNull) {
            return Double.valueOf(val);
        } else {
            throw new IllegalArgumentException("No value found for key " + key);
        }
    }
    
    public boolean getShiftCorrection(boolean defaultValue) {
        return Boolean.valueOf(props.getProperty("SHIFT_CORRECTION", defaultValue + ""));
    }

    public int getEdgeLearningType(int defaultValue, ModelType type) {
        String parameterName = "EDGE_LEARNING_TYPE";
        if (type == null) {
            // std value already set
        } else if (type == ModelType.FOREGROUND) {
            parameterName = "EDGE_LEARNING_TYPE_FG";
        } else if (type == ModelType.BACKGROUND) {
            parameterName = "EDGE_LEARNING_TYPE_BG";
        }
        return Integer.valueOf(props.getProperty(parameterName, defaultValue + ""));
    }

    public int getMotifLength(int defaultValue) {
        return Integer.valueOf(props.getProperty("MOTIF_LENGTH", defaultValue + ""));
    }

    public byte getMotifOrder(int defaultValue) {
        return Byte.valueOf(props.getProperty("MOTIF_ORDER", defaultValue + ""));
    }

    public byte getBackgroundOrder(int defaultValue) {
        return Byte.valueOf(props.getProperty("BACKGROUND_ORDER", defaultValue + ""));
    }

    public int getEMRestarts(int defaultValue) {
        return Integer.valueOf(props.getProperty("EM_RESTARTS", defaultValue + ""));
    }

    public double getProbThreshForWeights(double defaultValue) {
        return Double.valueOf(props.getProperty("PROB_THRESH_FOR_WEIGHTS", defaultValue + ""));
    }

    public double getMotifQualityThreshForWeights(double defaultValue) {
        return Double.valueOf(props.getProperty("MOTIF_QUALITY_THRESHOLD", defaultValue + ""));
    }

    /**
     * @param defaultValue
     * @param type
     * @return
     */
    public int learnFilterParameter(int defaultValue, ModelType type) {
        String parameterName = getParameterName("ASSUME_FILTERING", type);
        return Integer.valueOf(props.getProperty(parameterName, defaultValue + ""));
    }

    /**
     * 
     * @param defaultValue
     * @param type
     * @return
     */
    public int getNewickString(String defaultValue, ModelType type) {
        String parameterName = getParameterName("NEWICK", type);
        return Integer.valueOf(props.getProperty(parameterName, defaultValue + ""));
    }

    /**
     * A parameter can be set for background and/or foreground. If a parameter was set for background "_BG" is appended.
     * If a parameter was set for foreground then "_FG" is appended to the name.
     */
    private String getParameterName(String base, ModelType type) {
        if (type == null) {
            // std value already set
        } else if (type == ModelType.FOREGROUND) {
            base += ModelType.FOREGROUND.parameterSuffix;
        } else if (type == ModelType.BACKGROUND) {
            base += ModelType.BACKGROUND.parameterSuffix;
        }
        return base;
    }
}
