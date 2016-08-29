package tools;

import java.util.Properties;

import de.jstacs.data.AlphabetContainer;
import de.jstacs.data.Sample;
import de.jstacs.data.alphabets.UnobservableDNAAlphabet;
import evolution.EvolModel;
import io.SampleUtil;
import models.PhyloBackground;
import training.TrainingUtil;
import util.Config;

/**
 * This class is used to determine the branch lengths of a phylogenetic tree on a dataset using an phylogenetic
 * inhomogenious markov model of order 0.
 * 
 * @author Chaos
 * 
 */
public class LearnTree {
    public static void main(String[] args) throws Exception {
		EvolModel.MODEL_CLASS = "FS81alpha";
        Properties props = new Properties();
        Config.parseProperties(props, args, false);
        Config.replaceInternalPlaceholder(props);
        String newick = props.getProperty("model.newick");
        
        AlphabetContainer alphabet = new AlphabetContainer(new UnobservableDNAAlphabet());

        Sample[] samples = SampleUtil.getDataSets(props, alphabet);
        Sample union = Sample.union(samples[0],samples[1],samples[2],samples[3]);
        
        PhyloBackground.ENABLE_CACHING = true;
        PhyloBackground.EDGE_LEARNING_ONE_LENGTH = false;

		String learnedTree = TrainingUtil.learnTree(union, newick);

        System.out.println(learnedTree);
    }
}
