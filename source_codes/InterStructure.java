import java.util.Vector;
import cern.colt.list.DoubleArrayList;
import cern.jet.random.engine.MersenneTwister;
import cern.jet.random.engine.RandomEngine;
import cern.jet.random.sampling.RandomSampler;

public class InterStructure {

	int interN;	// # of nodes of an intercellular network 
	int interNPow; // 2^interN
	int interLink; // total # of links of an intercellular network
	boolean interAdjMatrix[][]; // adjacency matrix of an intercellular network  

	public void initialize(){
		
		// inter-network: node & link
		interN = 9; // # of cells		
		interNPow = (int) (Math.pow(interN,2)); 
		interLink = (int) (Math.random() * interNPow) + 1; // dynamic cellular topology: choose a random number between 1 and 81
		interAdjMatrix = new boolean[interN][interN]; 
	}

	public void makeStructure(){

		double[] nodeForLink = new double[interNPow];
		DoubleArrayList setForLink = new DoubleArrayList(nodeForLink);
		DoubleArrayList choiceForLink = new DoubleArrayList();
		for(int i = 0; i < interNPow; i++){
			nodeForLink[i] = i;
		}
		setForLink.shuffle();

		int Ls= setForLink.size();

		RandomEngine reS= new MersenneTwister(new java.util.Date());

		int link = interLink;
		long[] sampleForlink = new long[link];

		RandomSampler rsS = new RandomSampler(link, Ls, 0, reS);
		rsS.nextBlock(link, sampleForlink, 0);

		for (int j = 0; j < link; j++) {

			long num = sampleForlink[j];
			choiceForLink.add(setForLink.get((int) num));

		}

		for(int linkIdx = 0; linkIdx < interLink; linkIdx++){
			interAdjMatrix[(int) (choiceForLink.get(linkIdx) / interN)][(int) (choiceForLink.get(linkIdx) % interN)] = true;
		}
	}

	public void interStructure(){
		initialize();
		makeStructure();
	}
}
