package util;

/**
 * @author Martin Nettling
 * 
 */
public class IntPair implements Comparable<IntPair> {
    final public int val1;
    final public int val2;
    final public int allowedOverlapOfVal2;

    public IntPair(int v1, int v2) {
        this(v1, v2, 0);
    }

    public IntPair(int v1, int v2, int allowedOverlapOfVal2) {
        this.val1 = v1;
        this.val2 = v2;
        this.allowedOverlapOfVal2 = allowedOverlapOfVal2;
    }

    @Override
    public int compareTo(IntPair o) {
        if (val1 != o.val1) {
            return val1 - o.val1;
        } else {
            return val2 - o.val2;
        }
    }

    @Override
    public int hashCode() {
        return val1 + 31 * val2;
    }

    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        } else if (o instanceof IntPair &&
                ((IntPair) o).val1 == val1 && val2 == ((IntPair) o).val2) {
            return true;
        }
        return false;
    }

    public boolean softEquals(IntPair o) {
        if (o == this) {
            return true;
        } else if (((IntPair) o).val1 == val1 &&
                val2 - allowedOverlapOfVal2 <= ((IntPair) o).val2 &&
                ((IntPair) o).val2 <= val2 + allowedOverlapOfVal2) {
            return true;
        }
        return false;
    }

    @Override
    public String toString() {
        return val1 + "\t" + val2;
    }

    @Override
    public IntPair clone() {
        return new IntPair(val1, val2);
    }
}