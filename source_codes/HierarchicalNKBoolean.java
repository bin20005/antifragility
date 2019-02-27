import java.util.Arrays;
import java.util.Vector;

public class HierarchicalNKBoolean {

	InterStructure is = new InterStructure();
	NKBoolean nk = new NKBoolean();
	NKBoolean intraNet = null;

	int interN = 0;
	int interLink = 0;
	int hNBit = 0;	// intraNet n * interNet N
	int hPow = 0;	// 2^(n*N)

	static int hPowLimit = 10;// # of random initial states

	// Parameters on Antifragility----------------------------------	
	int simulTime_T = 400;	// simulation time T
	int perSize_X = 1;	// perturbed node size X
	int perPeriod_O = 1;	// perturbation frequency O		
	Vector<Double> delSigma = new Vector<Double>();	// vector to store mean of complexity differences
	Vector<Double> initialComplexityStorage = new Vector<Double>();	// vector of average complexity for all nodes before adding perturbations
	Vector<Double> finalComplexityStorage = new Vector<Double>();	// vector of average complexity for all nodes after adding perturbations 
	Vector<Double> fragility =  new Vector<Double>();	//vector of the value of fragility (-1<= fragility <=1). initial value = -2
	//--------------------------------------------------------------

	boolean[][] interAdjMatrix = null;  

	int[][] hNode = null;
	int[][] hNewNode = null;

	int n = 0;	
	int k = 0;
	int nPow = 0;

	Vector<Integer> communicatingNodeIDXVector = null;

	int[][] newDecNode = null;

	Vector<Integer> attractor = null;
	Vector<Integer> attrLeng = null; 
	int basinAttrSize[] = null;

	Vector<Integer>[] inputNodes = null;

	HierarchicalNKBoolean(HierarchicalNKBoolean Copier){

		this.interN = Copier.interN;
		this.interLink = Copier.interLink;
		this.hNBit = Copier.hNBit;
		this.hPow = Copier.hPow;

		this.n = Copier.n;	
		this.k = Copier.k;
		this.nPow = Copier.nPow;

		intVectorCopier(this.communicatingNodeIDXVector, Copier.communicatingNodeIDXVector);

		this.newDecNode = new int[nPow][1];

		this.attractor = new Vector<Integer>();
		this.attrLeng = new Vector<Integer>(); 
		this.basinAttrSize =  new int[Copier.attrLeng.size()];

		int2DArrayCopier(this.newDecNode, Copier.newDecNode, nPow, 1);
		intVectorCopier(this.attractor, Copier.attractor);
		intVectorCopier(this.attrLeng, Copier.attrLeng);
		int1DArrayCopier(this.basinAttrSize, Copier.basinAttrSize, attrLeng.size());

		this.interAdjMatrix = new boolean[interN][interN];  
		boolean2DArrayCopier(this.interAdjMatrix, Copier.interAdjMatrix, interN, interN);

		this.intraNet = Copier.intraNet;

	}

	HierarchicalNKBoolean(){
		initialize();
	}

	public void initialize(){

		nk.individualNetwork();
		intraNet = new NKBoolean(nk);

		is.interStructure();

		interN = is.interN;
		interLink = is.interLink;
		hNBit = 0;
		hPow = 0;

		interAdjMatrix = is.interAdjMatrix;  
		hNode = null;
		hNewNode = null;

		n = nk.n;
		k = nk.k;
		nPow = nk.nPow;

		communicatingNodeIDXVector = nk.communicatingNodeIDXVector;

		newDecNode = nk.newDecNode;
		attractor = nk.attractor;
		attrLeng = nk.attrLeng;
		basinAttrSize = nk.basinAttrSize;

		// to investigate inputs per node in an intercellular network
		inputNodes = new Vector[interN];

		for(int i = 0; i < interN; i++){
			inputNodes[i] = new Vector();
			for(int j = 0; j < interN; j++){
				if( interAdjMatrix[i][j] == true){
					inputNodes[i].add(j);
				}
			}
		}
	}

	public void detHierarchicalSTD(){

		hNBit = n * interN;

		// 2D array to store random initial states 
		int[][] randomInitial = new int[hPowLimit][1];

		// initialization for randomInitial array
		for(int ii = 0; ii < hPowLimit; ii++){
			randomInitial[ii] = null;
		}	

		// choosing random initial states (sampling without replacement)
		for(int i = 0; i < hPowLimit; i++){

			int[] tempRandom = new int[hNBit];

			for(int j = 0; j < hNBit; j++){
				tempRandom[j] = (int)(Math.random()*2);
			}		

			if(i == 0){
				randomInitial[i] = tempRandom;
			}
			else{

				int countCriterion = 0;
				int z = 0;
				while(true){
					if(randomInitial[z] != null){
						countCriterion++;
						z++;
					}
					else{
						break;
					}
				}

				int countEqual = 0;
				for(int w = 0; w < countCriterion; w++){
					if(!Arrays.equals(randomInitial[w], tempRandom)){
						countEqual++;
					}
				}

				if(countEqual == countCriterion){
					randomInitial[i] = tempRandom;
				}
				else{
					i--;
				}	
			}			
		}

		// tracking state transitions until T
		int[][] configurationTrack = new int[(simulTime_T*2)+1][1]; 
		int[][] configurationTrack2 = new int[(simulTime_T*2)+1][1]; 

		// before adding perturbations
		int[] tempConfigBinary = new int[hNBit];	
		int[] tempNewConfigBinary = new int[hNBit];

		// after adding perturbations
		int[] tempConfigBinary2 = new int[hNBit];	
		int[] tempNewConfigBinary2 = new int[hNBit];
		double tempInitialComplexity;
		double tempFinalComplexity;

		for(int conIdx = 0; conIdx < hPowLimit; conIdx++){

			// initial & final complexity initialization
			tempInitialComplexity = -2;
			tempFinalComplexity = -2;

			// initialization of configurationTrack array
			for(int ii = 0; ii < simulTime_T*2; ii++){
				configurationTrack[ii] = null;
				configurationTrack2[ii] = null;	
			}

			configurationTrack[0] = randomInitial[conIdx];
			configurationTrack2[0] = randomInitial[conIdx];	// same random initial states

			for(int init1 = 0; init1 < hNBit; init1++){	// initialization for tempConfigBinary array
				tempConfigBinary[init1] = randomInitial[conIdx][init1];
				tempConfigBinary2[init1] = randomInitial[conIdx][init1];
			}			

			for(int init2 = 0; init2 < hNBit; init2++){	// initialization for tempNewConfigBinary array
				tempNewConfigBinary[init2] = -1;
				tempNewConfigBinary2[init2] = -1;
			}

			int time = 0;

			while(time < simulTime_T*2){

				for(int hNodeIdx = 0; hNodeIdx < interN; hNodeIdx++){

					if(inputNodes[hNodeIdx] != null){	// if there are inputs

						// (i) update of communicating nodes
						int communicatingNodeIDX = -1;
						for(int comm = 0; comm < communicatingNodeIDXVector.size(); comm++){

							communicatingNodeIDX = communicatingNodeIDXVector.get(comm);

							Vector<Integer>[] compareNodeVal =  new Vector[inputNodes[hNodeIdx].size()];
							Vector<Integer>[] compareNodeVal2 =  new Vector[inputNodes[hNodeIdx].size()];
							for(int initIdx = 0; initIdx < inputNodes[hNodeIdx].size(); initIdx++){
								compareNodeVal[initIdx] = new Vector();
								compareNodeVal2[initIdx] = new Vector();
							}

							for(int itrInputIdx = 0; itrInputIdx < inputNodes[hNodeIdx].size(); itrInputIdx++){
								int tempInputNode = inputNodes[hNodeIdx].get(itrInputIdx);

								int tempBinVal = tempConfigBinary[(n * tempInputNode) + communicatingNodeIDX];
								compareNodeVal[itrInputIdx].add(tempBinVal);

								int tempBinVal2 = tempConfigBinary2[(n * tempInputNode) + communicatingNodeIDX];
								compareNodeVal2[itrInputIdx].add(tempBinVal2);
							}
 
							// original RBN
							int oneCount = 0;
							for(int compareIdx = 0; compareIdx < inputNodes[hNodeIdx].size(); compareIdx++){
								if(compareNodeVal[compareIdx].get(0) == 1){
									oneCount++;
								}
							}

							if(oneCount > 0){
								tempNewConfigBinary[(n*hNodeIdx) + communicatingNodeIDX] = 1;
							}
							else{
								tempNewConfigBinary[(n*hNodeIdx) + communicatingNodeIDX] = 0;
							}

							//perturbed RBN
							int oneCount2 = 0;
							for(int compareIdx2 = 0; compareIdx2 < inputNodes[hNodeIdx].size(); compareIdx2++){
								if(compareNodeVal2[compareIdx2].get(0) == 1){
									oneCount2++;
								}
							}

							if(oneCount2 > 0){
								tempNewConfigBinary2[(n*hNodeIdx) + communicatingNodeIDX] = 1;
							}
							else{
								tempNewConfigBinary2[(n*hNodeIdx) + communicatingNodeIDX] = 0;
							}
						}


						// (ii) update of the other nodes
						// original RBN
						int tempDecVal = 0;
						int p = n - 1;
						for(int bitIdx = 0; bitIdx < n; bitIdx++){
							tempDecVal = tempDecVal + tempConfigBinary[(n*hNodeIdx)+bitIdx] * (int) Math.pow(2,p);
							p--;
						}

						tempDecVal = newDecNode[tempDecVal][0];

						for(int biTransIdx = n - 1; biTransIdx >= 0; biTransIdx--){
							if(!communicatingNodeIDXVector.contains(biTransIdx)){
								tempNewConfigBinary[(n*hNodeIdx)+biTransIdx] = tempDecVal % 2;
							}
							tempDecVal = tempDecVal / 2;
						}

						// perturbed RBN
						int tempDecVal2 = 0;
						int p2 = n - 1;
						for(int bitIdx2 = 0; bitIdx2 < n; bitIdx2++){
							tempDecVal2 = tempDecVal2 + tempConfigBinary2[(n*hNodeIdx)+bitIdx2] * (int) Math.pow(2,p2);
							p2--;
						}

						tempDecVal2 = newDecNode[tempDecVal2][0];

						for(int biTransIdx2 = n - 1; biTransIdx2 >= 0; biTransIdx2--){
							if(!communicatingNodeIDXVector.contains(biTransIdx2)){
								tempNewConfigBinary2[(n*hNodeIdx)+biTransIdx2] = tempDecVal2 % 2;
							}
							tempDecVal2 = tempDecVal2 / 2;
						}
					}
					else{	// if there is no input

						// original RBN
						int tempDecVal = 0;
						int p = n - 1;
						for(int bitIdx = 0; bitIdx < n; bitIdx++){
							tempDecVal = tempDecVal + tempConfigBinary[(n*hNodeIdx)+bitIdx] * (int) Math.pow(2,p);
							p--;
						}

						tempDecVal = newDecNode[tempDecVal][0];

						for(int biTransIdx = n - 1; biTransIdx >= 0; biTransIdx--){
							tempNewConfigBinary[(n*hNodeIdx)+biTransIdx] = tempDecVal % 2;
							tempDecVal = tempDecVal / 2;
						}

						// perturbed RBN
						int tempDecVal2 = 0;
						int p2 = n - 1;
						for(int bitIdx2 = 0; bitIdx2 < n; bitIdx2++){
							tempDecVal2 = tempDecVal2 + tempConfigBinary2[(n*hNodeIdx)+bitIdx2] * (int) Math.pow(2,p2);
							p2--;
						}

						tempDecVal2 = newDecNode[tempDecVal2][0];

						for(int biTransIdx2 = n - 1; biTransIdx2 >= 0; biTransIdx2--){
							tempNewConfigBinary2[(n*hNodeIdx)+biTransIdx2] = tempDecVal2 % 2;
							tempDecVal2 = tempDecVal2 / 2;
						}
					}
				}

				int[] temp = new int[hNBit];
				int[] temp2 = new int[hNBit];

				int1DArrayCopier(temp,tempNewConfigBinary, hNBit);
				int1DArrayCopier(temp2,tempNewConfigBinary2, hNBit);
				time++;

				configurationTrack[time] = temp;	
				if(time%perPeriod_O == 0){
					// choosing perturbed nodes
					Vector<Integer> randNodeIDX4Per = new Vector<Integer>();
					for(int p = 0; p < perSize_X; p++){
						int tempRandIDX = (int)(Math.random()*hNBit);
						randNodeIDX4Per.addElement(tempRandIDX);
					}

					for(int f = 0; f < randNodeIDX4Per.size(); f++){
						int idx = randNodeIDX4Per.get(f);
						// perturbation: flip all the states of chosen nodes
						if(temp2[idx] == 0){
							temp2[idx] = 1;
						}
						else{
							temp2[idx] = 0;
						}
					}

					configurationTrack2[time] = temp2;

				}
				else{
					configurationTrack2[time] = temp2;	
				}
				int1DArrayCopier(tempConfigBinary,temp, hNBit);
				int1DArrayCopier(tempConfigBinary2,temp2, hNBit);

			}

			// measuring complexity based on Shannon's entropy 
			double totalE = 0;
			double totalE2 = 0;
			double initialComplexity = 0;
			double finalComplexity = 0;
			for(int nIDX = 0; nIDX < hNBit; nIDX++){
				int p0 = 0;
				int p1 = 0;
				int p0_2 = 0;
				int p1_2 = 0;

				for(int c = simulTime_T+1; c < (simulTime_T*2)+1; c++){
					if(configurationTrack[c][nIDX] == 0){
						p0++;
					}
					else{
						p1++;
					}

					if(configurationTrack2[c][nIDX] == 0){
						p0_2++;
					}
					else{
						p1_2++;
					}

				}	

				double E = -1;
				if(p0 == 0 | p1 == 0){
					E = 0;
				}
				else{
					double prob = (double) p0/simulTime_T;
					E = -(prob*log(2, prob) + (1-prob)*log(2, (1-prob))); 
				}

				double E2= - 1;
				if(p0_2 == 0 | p1_2 == 0){
					E2 = 0; 
				}
				else{
					double prob2 = (double) p0_2/simulTime_T;
					E2 = -(prob2*log(2, prob2) + (1-prob2)*log(2, (1-prob2)));
				}
				totalE += E;
				totalE2 += E2;
			}
			initialComplexity = 4*((double)totalE/hNBit)*(1-((double)totalE/hNBit));
			finalComplexity = 4*((double)totalE2/hNBit)*(1-((double)totalE2/hNBit));
			
			initialComplexityStorage.addElement(initialComplexity);
			finalComplexityStorage.addElement(finalComplexity);
			delSigma.addElement(finalComplexity - initialComplexity);
		}
	}

	public void calFragility(){
		
		double degreeOfPer = (perSize_X*((double) simulTime_T/perPeriod_O))/((double) hNBit*simulTime_T);
		double tempFragility = -2;
		for(int cv = 0; cv < delSigma.size(); cv++){
			tempFragility = -delSigma.get(cv) * degreeOfPer;
			fragility.addElement(tempFragility);
		}
	}

	public void hierarchicalNetwork(){
		detHierarchicalSTD();
		calFragility();
	}

	double log(double base, double num)
	{
		return Math.log(num) / Math.log(base);
	}

	public void intVectorCopier(Vector<Integer> targetVector, Vector<Integer> sourceVector)
	{
		Integer x;
		for(int i=0; i < sourceVector.size(); i++)
		{
			x = new Integer(sourceVector.get(i).intValue());
			targetVector.add( x );
		}
	}

	public void longVectorCopier(Vector<Long> targetVector, Vector<Long> sourceVector)
	{
		Long x;
		for(int i=0; i < sourceVector.size(); i++)
		{
			x = new Long(sourceVector.get(i).longValue());
			targetVector.add( x );
		}
	}

	public void boolean2DArrayCopier(boolean targetArray[][], boolean sourceArray[][], int x, int y)
	{
		for(int i=0; i < x; i++)
		{
			for(int j=0; j < y; j++)
			{
				targetArray[i][j] = sourceArray[i][j];
			}
		}
	}

	public void int1DArrayCopier(int targetArray[], int sourceArray[], int x)
	{
		for(int i=0; i < x; i++)
		{
			targetArray[i] = sourceArray[i];
		}
	}

	public void long1DArrayCopier(long targetArray[], long sourceArray[], int x)
	{
		for(int i=0; i < x; i++)
		{
			targetArray[i] = sourceArray[i];
		}
	}

	public void double1DArrayCopier(double targetArray[], double sourceArray[], int x)
	{
		for(int i=0; i < x; i++)
		{
			targetArray[i] = sourceArray[i];
		}
	}

	public void int2DArrayCopier(int targetArray[][], int sourceArray[][], int x,int y)
	{
		for(int i=0; i < x; i++)
		{
			for(int j=0; j < y; j++)
			{
				targetArray[i][j] = sourceArray[i][j];
			}
		}
	}

	public void long2DArrayCopier(long targetArray[][], long sourceArray[][], int x,int y)
	{
		for(int i=0; i < x; i++)
		{
			for(int j=0; j < y; j++)
			{
				targetArray[i][j] = sourceArray[i][j];
			}
		}
	}

	public void intVectorArrayCopier(Vector<Integer>[] targetVectorArray, Vector<Integer>[] sourceVectorArray)
	{
		for(int i = 0; i < interN; i++){
			targetVectorArray[i] = new Vector<Integer>();
			for(int j = 0; j < sourceVectorArray[i].size(); j++){
				targetVectorArray[i].add(sourceVectorArray[i].get(j));
			}

		}
	}

	public void intVectorArrayCopier2(Vector<Integer>[] targetVectorArray, Vector<Integer>[] sourceVectorArray, int x)
	{
		for(int i = 0; i < x; i++){
			targetVectorArray[i] = new Vector<Integer>();
			for(int j = 0; j < sourceVectorArray[i].size(); j++){
				targetVectorArray[i].add(sourceVectorArray[i].get(j));
			}

		}
	}

	public void longVectorArrayCopier(Vector<Long>[] targetVectorArray, Vector<Long>[] sourceVectorArray, int x)
	{
		for(int i = 0; i < x; i++){
			targetVectorArray[i] = new Vector<Long>();
			for(int j = 0; j < sourceVectorArray[i].size(); j++){
				targetVectorArray[i].add(sourceVectorArray[i].get(j));
			}

		}
	}

	public void doubleVectorArrayCopier(Vector<Double>[] targetVectorArray, Vector<Double>[] sourceVectorArray, int x)
	{
		for(int i = 0; i < x; i++){
			targetVectorArray[i] = new Vector<Double>();
			for(int j = 0; j < sourceVectorArray[i].size(); j++){
				targetVectorArray[i].add(sourceVectorArray[i].get(j));
			}

		}
	}

}
