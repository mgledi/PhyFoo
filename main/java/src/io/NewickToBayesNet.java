package io;

import java.io.IOException;
import java.io.Reader;
import java.io.StreamTokenizer;
import java.io.StringReader;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import bayesNet.BayesNet;
import bayesNet.BayesNetNode;

/**
 * This class contains static methods to create a {@link BayesNet} from a String in newick format
 * 
 * @author Martin Nettling
 *
 */
public class NewickToBayesNet {

    private static String prefix;

    public static BayesNet getTree(String r, String name) {
        if (!r.endsWith(";")) {
            r += ";";
        }
        prefix = name;
        StringReader sr = new StringReader(r);
        return getTree(sr, name);
    }

    public static BayesNet getTree(String r) {
        prefix = "default";
        return getTree(r, "default");
    }

    public static BayesNet getTree(StringReader r, String name) {
        NewickStreamTokenizer st = new NewickStreamTokenizer(r);
        BayesNet t = new BayesNet(name);
        if (t.getRoot() == null) {
            // t.setRoot(new BayesNetNode(t.getName(),t, BayesNetNode.PHYLO_ROOT_NODE));
            t.setRoot(new BayesNetNode(t.getName(), t));
            t.getRoot().props.setPhyloRoot();
        }

        StringBuffer label = new StringBuffer();
        StringBuffer nhx = new StringBuffer();

        BayesNetNode n = t.getRoot();
        BayesNetNode c = null;

        double d = -1.0;
        boolean read = true;
        boolean distance = false;
        boolean addNhx = false;
        boolean escaped = false;
        boolean somethingRead = false;
        boolean finished = false;

        while (read) {
            try {
                int next = st.nextToken();
                if (addNhx && next != ']') {
                    if (st.sval == null) {
                        nhx.append((char) next);
                    }
                    if (st.sval != null) {
                        nhx.append(st.sval);
                    }
                    continue;
                }

                if (escaped) {
                    if (st.sval != null) {
                        label.append(st.sval);
                    } else {
                        label.append(Character.toString((char) next));
                    }
                    escaped = false;
                    continue;
                }

                switch (next) {
                // ############ es wird ein Wort gelesen ( LABEL oder Distanz oder Kommentar )
                case StreamTokenizer.TT_WORD:
                    if (distance) {
                        d = Double.parseDouble(st.sval);
                        distance = false;
                    } else {
                        label.append(st.sval);
                    }
                    break;
                // ############ Stream ist am Ende
                case StreamTokenizer.TT_EOF: {
                    read = false;
                    break;
                }
                    // ############ einzelF�lle des StringTokenizers richtig behandeln
                default:
                    // solange kein Semicolon, solange wurde der String noch nicht abgeschlossen
                    if (st.ttype != ';') {
                        somethingRead = true;
                    }

                    switch (st.ttype) {
                    // neuer Subtree wird angelegt
                    case '(':
                        // c = new BayesNetNode(t, BayesNetNode.PHYLO_INNER_NODE); // neuen Knoten anlegen, soll Kind
                        // von aktuellem Knoten werden
                        c = new BayesNetNode(t); // neuen Knoten anlegen, soll Kind von aktuellem Knoten werden
                        if (t.getRoot() == null) {
                            // t.setRoot(c);
                        }

                        n.addChild(c);
                        c.setName(prefix + "iN" + c.nodeNumber);
                        c.setRealName(prefix);
                        n = c;
                        break;
                    // Subtree wird geschlossen
                    case ')':
                        setLabel(n, label, nhx);
                        if (d != -1.0) {
                            n.setDistanceToParent(d);
                            d = -1.0;
                        }
                        label = new StringBuffer(); // reinitialisiere StringBuffer
                        nhx = new StringBuffer(); // reinitialisiere StringBuffer
                        n = n.getParent(0);
                        break;
                    case ',':
                        setLabel(n, label, nhx);
                        if (d != -1.0) {
                            n.setDistanceToParent(d);
                            d = -1.0;
                        }

                        label = new StringBuffer(); // reinitialisiere StringBuffer
                        nhx = new StringBuffer(); // reinitialisiere StringBuffer
                        // c = new BayesNetNode(t, BayesNetNode.PHYLO_INNER_NODE);
                        c = new BayesNetNode(t);
                        c.setName(prefix + "_iN" + c.nodeNumber);
                        c.setRealName(prefix);
                        n.getParent(0).addChild(c);
                        n = c;
                        break;
                    // ####### lese die Distanz ein
                    case ':':
                        distance = true;
                        break;
                    case '[':
                        int peek = st.nextToken();
                        if (peek == ']') {
                            addNhx = true;
                        } else {
                            System.out.println("no ! not nhx");
                            label.append("[");
                        }
                        st.pushBack();

                        break;
                    case ']':
                        if (addNhx) {
                            addNhx = false;
                        } else {
                            label.append("]");
                        }
                        break;
                    // ####### Ende des Strings schreibe letzten Kntoen
                    case ';':
                        // n.setName(label.toString());
                        if (d != -1.0) {
                            // sollte Baum nur aus wurzel Bestehen
                            if (!n.isRoot()) {
                                // //TODO: Distance hinzuf�gen
                                // n.setDistanceToParent(d);
                                d = -1.0;
                            }
                        }
                        label = new StringBuffer(); // reinitialisiere StringBuffer
                        nhx = new StringBuffer(); // reinitialisiere StringBuffer
                        read = false; // Stream zu Ende gelesen
                        finished = true; // und abgeschlossen
                        break;
                    case '\'':
                    case '"':
                        label.append(st.sval);
                        break;
                    case '\\':
                        escaped = true;
                        break;
                    }
                    break;
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        if (!finished && somethingRead) {
            System.err.println("Error reading the newickstring to a tree.");
        }
        for(BayesNetNode node : t.nodes) {
            if(node.isLeaf()) node.props.setPhyloLeaf();
        }
        return t;
    }

    private static void setLabel(BayesNetNode n, StringBuffer label, StringBuffer nhxB) {
        if (nhxB.length() > 0) {
            String nhx = nhxB.substring(6);
            String[] tags = nhx.split(":");
            HashMap<String, String> extLabel = new HashMap<String, String>();
            for (int k = 0; k < tags.length; k++) {
                String[] nameValue = tags[k].split("=");
                extLabel.put(nameValue[0], nameValue[1]);
            }
            // TODO: Extra Taxoninformationen mit im Baum speichern
        }
        if (label.length() > 0) {
            String l = prefix + "_" + label.toString();
            l = l.trim();
            n.setName(l);
            n.setRealName(label.toString());
        }
    }

    public static String getStringFromTree(BayesNet t) {
        return getStringFromTree(t, true);
    }

    public static String getStringFromTree(BayesNet t, boolean addSemicolon) {
        StringBuffer treeString = new StringBuffer();
        BayesNetNode root = t.getRoot();
        appendNode(root, treeString);
        return treeString.toString() + (addSemicolon ? ";" : "");
    }

    /**
     * This is an internal method to parse the given String buffer. The method converts the given node to newick string
     * and appends the result to the given StringBuffer.
     * 
     * @param node
     * @param treeString
     */
    private static void appendNode(BayesNetNode node, StringBuffer s) {
        if (node == null)
            return;

        if (!node.isLeaf()) {
            s.append("(");
            boolean added = false;
            for (BayesNetNode n : node.children) {
                appendNode(n, s);
                s.append(',');
                added = true;
            }
            if (added)
                s.deleteCharAt(s.length() - 1);
            s.append(")");
        }
        if (node.getName() != null) {
            /*
             * Escape characters in the label
             * 
             * Escape all : as \: Escape all " " (space) as "\ "
             */

            String label = node.getName();
            if (label.indexOf(":") >= 0) {
                Pattern p = Pattern.compile(":");
                Matcher matcher = p.matcher(label);
                ;
                label = matcher.replaceAll("\\\\$0");
            }

            Pattern spacePattern = Pattern.compile("\\s");
            Matcher matcher = spacePattern.matcher(label);
            ;
            label = matcher.replaceAll("\\\\$0");

            s.append(label);
        }

        // TODO: Distance hinzuf�gen
        // if(node.getDistanceToParent() != 1.0 && node.getDistanceToParent() >= 0)
        // {
        // s.append(":" + node.getDistanceToParent());
        // }
        // TODO: extra Annotationen mitschreiben
        // writeNHXExtensiont(node, s);
    }
}

final class NewickStreamTokenizer extends StreamTokenizer {
    public NewickStreamTokenizer(Reader r) {
        super(r);
        resetSyntax();
        wordChars(0, 255);
        whitespaceChars(0, '\n');

        ordinaryChar(';');
        ordinaryChar(',');
        ordinaryChar(')');
        ordinaryChar('(');
        ordinaryChar('[');
        ordinaryChar(']');
        ordinaryChar(':');
        ordinaryChar('\\');
        // commentChar('/');
        quoteChar('"');
        quoteChar('\'');
    }
}
