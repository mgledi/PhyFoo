package util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Map.Entry;
import java.util.Arrays;
import java.util.Properties;
import java.util.logging.Logger;

/**
 * This class helps to configure PhyloBayes.
 * 
 * @author Martin Nettling
 * 
 */
public class Config {
    static Logger LOG = Logger.getLogger(Config.class.getName());

    /**
     * The base directory. All other used directories are assumed to be here, if not given absolute.
     */
    public static String BASE_DIR = "./";

    static {
        if (new File("C:/Users/Martin/Dropbox/promotion/").exists()) {
            BASE_DIR = "E:/Dropbox/promotion/";
        } else if (new File("E:/Dropbox/promotion/").exists()) {
            BASE_DIR = "E:/Dropbox/promotion/";
        } else if (new File("/data/m.gleditzsch/Dropbox/promotion/").exists()) {
            BASE_DIR = "/data/m.gleditzsch/Dropbox/promotion/";
        } else if (new File("/home/gleditzsch/").exists()) {
            BASE_DIR = "/home/gleditzsch/";
        }
    }

    /**
     * @param filename
     *            the path to the file, that contains the properties
     * @param step
     * @return the retrieved properties
     * @throws FileNotFoundException
     * @throws IOException
     */
    public static Properties getProperties(String filename, int step) throws FileNotFoundException, IOException {
        Properties p = new Properties();
        p.load(new FileInputStream(new File(filename)));
        for (Entry<Object, Object> entry : p.entrySet()) {
            String key = (String) entry.getKey();
            String value = (String) entry.getValue();
            value = value.replace("%step%", "" + step);
            p.setProperty(key, value);
        }
        return p;
    }

    /**
     * @param props
     * @param args
     */
    public static void parseProperties(Properties props, String[] args) {
        parseProperties(props, args, true);
    }

    /**
     * Parses the properties given in the array args and copies them to props. If keyMustExist is TRUE than only
     * properties that already exists will be respected.
     * 
     * @param props
     * @param args
     * @param keyMustExist
     */
    public static void parseProperties(Properties props, String[] args, boolean keyMustExist) {
        for (int i = 0; i < args.length; i++) {
            if (args[i].startsWith("--")) {
                String tmp = args[i].substring(2);
                String key = tmp.substring(0, tmp.indexOf('='));
                String value = tmp.substring(tmp.indexOf('=') + 1);

                if (props.containsKey(key) || !props.containsKey(key) && !keyMustExist) {
                    props.setProperty(key, value);
                    if (keyMustExist && !props.containsKey(key)) {
                        LOG.warning("Can't find a property with the key: " + key + " However, property was set");
                    }
                } else {
                    throw new IllegalArgumentException("Can't find a property with the key: " + key
                            + ". The key must be of " + props.keySet());
                }
            }
        }
    }

    /**
     * Copies the properties in from to the properties in to. Properties in to with the same name will be overwritten.
     * 
     * @param from
     * @param to
     */
    public static void copyProperties(Properties from, Properties to) {
        for (Object key : from.keySet()) {
            to.setProperty((String) key, (String) from.get(key));
        }
    }

    /**
     * Tries to access the property with the given key in the given props. If allowNull is TRUE this method will return
     * null if the key does not exist or throws an {@link IllegalArgumentException}.
     * 
     * @param props
     *            the properties of interest
     * @param key
     *            the key to access
     * @param allowNull
     *            if true, the method is allowed to return null
     * @return the corresponding value or null
     */
    public static Value getProperty(Properties props, String key, boolean allowNull) {
        String val = props.getProperty(key);
        if (val != null || allowNull) {
            return new Value(val);
        } else {
            throw new IllegalArgumentException("No value found for key \"" + key + "\"");
        }
    }

    /**
     * Tries to access the property with the given key in the given props. If no value can be found for the given key,
     * the value for the fallbackKey is tried to access. If allowNull is TRUE and no value was found this method will
     * return
     * null or throws an {@link IllegalArgumentException} anyway.
     * 
     * @param props
     *            the properties of interest
     * @param allowNull
     *            if true, the method is allowed to return null
     * @param keys
     *            the keys to access. The first key thats value isn't null is used.
     * @return the corresponding value or null
     */
    public static Value getPropertyWithFallBack(Properties props, boolean allowNull, String... keys) {
        String val = null;
        int i = 0;
        while (val == null && i < keys.length) {
            val = props.getProperty(keys[i++]);
        }
        if (val != null || allowNull) {
            return new Value(val);
        } else {
            throw new IllegalArgumentException("No value found for key \"" + Arrays.toString(keys) + "\"");
        }
    }

    /**
     * Tries to access the property with the given key in the given props. If no value can be found for the given key,
     * the given defaultValue is returned.
     * 
     * @param props
     *            the properties of interest
     * @param key
     *            the key to access
     * @param defaultValue
     *            the defaultValue to set if no value was found for key
     * @return the corresponding value or the defaultValue
     */
    public static Value getProperty(Properties props, String key, String defaultValue) {
        return new Value(props.getProperty(key, defaultValue));
    }

    /**
     * Tries to access the property with the given key in the given props. If no value can be found for the given key,
     * the value for the fallbackKey is tried to access. If again no value was found the given defaultValue is returned.
     * 
     * @param props
     *            the properties of interest
     * @param defaultValue
     *            the defaultValue to set if no value was found for key or fallbackKey
     * @param keys
     *            the key to access
     * @return the corresponding value or the defaultValue
     */
    public static Value getPropertyWithFallBack(Properties props, String defaultValue, String... keys) {
        Value value = getPropertyWithFallBack(props, true, keys);
        if (value.val == null) {
            return new Value(defaultValue);
        }
        return value;
    }

    /**
     * This is a very simple method to replace placeholders with internal values. Example:
     * 
     * <pre>
     * basefolder = myfolder
     * storefolder = %basefolder%/subfolder
     * </pre>
     * 
     * This method will replace %basefolder% with myfolder.
     * 
     * @param props
     *            your properties
     */
    public static void replaceInternalPlaceholder(Properties props) {
        for (Object key : props.keySet()) {
            String skey = (String) key;
            for (Object key2 : props.keySet()) {
                String skey2 = (String) key2;
                String sval2 = props.getProperty(skey2);
                if (sval2.contains("%")) {
                    sval2 = sval2.replace("%" + skey + "%", props.getProperty(skey));
                    props.setProperty(skey2, sval2);
                }
            }
        }
    }

    /**
     * This class is a wrapper for values in config. It allows easy converting by fluent style pattern.
     * 
     * @author Martin Nettling
     * 
     */
    public static class Value {
        private final String val;

        /**
         * Instantiates a new Value with the given object.
         * 
         * @param val
         */
        public Value(String val) {
            this.val = val;
        }

        /**
         * Casts the internal value to a string
         * 
         * @return the value as String
         */
        public String asString() {
            return val;
        }

        /**
         * Casts the internal value to a boolean
         * 
         * @return the value as boolean
         */
        public Boolean asBoolean() {
            return Boolean.valueOf(val);
        }

        /**
         * Casts the internal value to a float
         * 
         * @return the value as float
         */
        public float asFloat() {
            return Float.valueOf(val);
        }

        /**
         * Casts the internal value to a double
         * 
         * @return the value as double
         */
        public double asDouble() {
            return Double.valueOf(val);
        }

        /**
         * Casts the internal value to an int
         * 
         * @return the value as int
         */
        public int asInt() {
            return Integer.valueOf(val);
        }

        /**
         * Casts the internal value to a byte
         * 
         * @return the value as byte
         */
        public byte asByte() {
            return Byte.valueOf(val);
        }
        
        @Override
        public String toString() {
            return String.valueOf(val);
        }
    }
}
