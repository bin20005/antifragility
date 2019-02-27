import java.util.Arrays;
import java.util.Vector;

public class NKBoolean {

	int n,k; 							// # of nodes, # of links per node
	int nPow; 							// 2^n 
	int kPow;							// 2^k
	int kPowMax;						// 2^(n+1)
	Vector<Integer>[] appliedNode;	    // inputs per node
	Vector<Integer> attractor; 			// attractor states
	Vector<Integer> attrLeng; 			// attractor length
	int basinAttr[][]; 					// attractor each configuration converges into
	boolean booleFunc[][]; 				// Boolean functions per node
//	Vector<Boolean>[] booleFunc;		// Boolean functions of CD4+ T-cell 
	boolean adjMatrix[][]; 				// adjacency matrix
	boolean nodeVal[][]; 				// state transition table at time T (binary)
	boolean newNodeVal[][];				// state transition table at time T+1 (binary)
	int decNode[][];					// state transition table at time T (decimal)
	int newDecNode[][]; 				// state transition table at time T+1 (decimal)

	int attrLengIn;						// attractor length index 
	int basinAttrSize[];				// attractor length + basin size per attractor
	int onlyBasinSize[];				// basin size per attractor

	static int numOfCommunicatingNodes = 6;	// # of communicating nodes
	Vector<Integer> communicatingNodeIDXVector;	// communicating node index

	static int powLimit = 10;// # of random initial states, powLimit of CD4+ T-cell = 1000
	
	// Parameters on Antifragility----------------------------------
	int nodeLimit = 20;
	int simulTime_T =  400; // simulation time T
	int perSize_X = 2;	// perturbed node size X
	int perPeriod_O = 1;	// perturbation frequency O			
	Vector<Double> delSigma = new Vector<Double>();	// vector to store mean of complexity differences 
	Vector<Double> initialComplexityStorage = new Vector<Double>();
	Vector<Double> finalComplexityStorage = new Vector<Double>();
	Vector<Double> fragility =  new Vector<Double>();	//vector to store the value of fragility (-1<= fragility <=1). initial value = -2
	//--------------------------------------------------------------
	
	NKBoolean(NKBoolean Copier){

		this.n = Copier.n;
		this.k = Copier.k;

		this.kPow = Copier.kPow;
		this.kPowMax = Copier. kPowMax;
		
		this.booleFunc = new boolean[kPow][n]; 
//		this.booleFunc =  new Vector[n]; // CD4+ T-cell
		this.adjMatrix = new boolean[n][n]; 

		this.appliedNode =  new Vector[n];	   
		
		this.communicatingNodeIDXVector = new Vector<Integer>(); 

		intVectorCopier(this.communicatingNodeIDXVector, Copier.communicatingNodeIDXVector);

		intVectorArrayCopier(this.appliedNode, Copier.appliedNode); 

		boolean2DArrayCopier(this.booleFunc, Copier.booleFunc, kPow, n);
//		booleanVectorArrayCopier(this.booleFunc, Copier.booleFunc); // CD4+ T-cell
		boolean2DArrayCopier(this.adjMatrix, Copier.adjMatrix, n, n);


		if(n <= nodeLimit){
			this.nPow = Copier.nPow; 

			this.attractor = new Vector<Integer>();
			this.attrLeng = new Vector<Integer>(); 

			this.basinAttr = new int[nPow][3];
			this.nodeVal = new boolean[nPow][n]; 
			this.newNodeVal = new boolean[nPow][n]; 
			this.decNode = new int[nPow][1]; 
			this.newDecNode = new int[nPow][1]; 

			this.basinAttrSize = new int[Copier.attrLeng.size()];
			this.onlyBasinSize = new int[Copier.attrLeng.size()];

			intVectorCopier(this.attractor, Copier.attractor);
			intVectorCopier(this.attrLeng, Copier.attrLeng);

			int1DArrayCopier(this.basinAttrSize, Copier.basinAttrSize, attrLeng.size());
			int1DArrayCopier(this.onlyBasinSize, Copier.onlyBasinSize, attrLeng.size());

			int2DArrayCopier(this.basinAttr, Copier.basinAttr, nPow, 3);
			int2DArrayCopier(this.decNode, Copier.decNode, nPow, 1);
			int2DArrayCopier(this.newDecNode, Copier.newDecNode, nPow, 1);

			boolean2DArrayCopier(this.nodeVal, Copier.nodeVal, nPow, n);
			boolean2DArrayCopier(this.newNodeVal, Copier.newNodeVal, nPow, n);

			this.attrLengIn = Copier.attrLengIn; 
		}
	}	

	NKBoolean(){	
		initialize();
	}

	public void initialize()	
	{	
		// intracellular network: node & link
		n = 18;
		k = 1; // when the network is CD4+ T-cell, it is commented out
		
		communicatingNodeIDXVector = new Vector<Integer>();
		// choosing communicating nodes randomly
		int temp = -1;
		while(true){
			if(communicatingNodeIDXVector.size() == numOfCommunicatingNodes){
				break;
			}

			temp = (int) (Math.random() * n);
			if(!communicatingNodeIDXVector.contains(temp)){
				communicatingNodeIDXVector.add(temp);
			}
		}

//		// CD4+ T-cell
//		communicatingNodeIDXVector = new Vector<Integer>();
//		communicatingNodeIDXVector.add(12);
//		communicatingNodeIDXVector.add(13);
//		communicatingNodeIDXVector.add(14);
//		communicatingNodeIDXVector.add(15);
//		communicatingNodeIDXVector.add(16);
//		communicatingNodeIDXVector.add(17);
		
		kPow = (int) Math.pow(2,k);
		kPowMax = (int) Math.pow(2, k+1);

		booleFunc = new boolean[kPow][n];
//		booleFunc = new Vector[n]; // CD4+ T-cell
		adjMatrix = new boolean[n][n]; 

		appliedNode = new Vector[n];	    


		if(n <= nodeLimit){
			nPow = (int) Math.pow(2,n); 

			nodeVal = new boolean[nPow][n]; 
			newNodeVal = new boolean[nPow][n]; 
			decNode = new int[nPow][1]; 
			newDecNode = new int[nPow][1]; 

			attractor = new Vector<Integer>(); 
			basinAttr = new int[nPow][3]; 
			attrLeng = new Vector<Integer>(); 

			basinAttrSize = null;
			onlyBasinSize = null;	
			attrLengIn = 0;
		}
	}

	public void makeStructure(){

		// adjacency matrix
		for(int xx = 0; xx < n; xx++){
			for(int yy = 0; yy < n; yy++){
				adjMatrix[xx][yy] = false;
			}
		}

		// random topology
		for(int x = 0;  x < n; x++){

			Vector <Integer> inputs = new Vector<Integer>(); 
			int randNum = -1;
			for(int numOfInputs = 0; numOfInputs < k; numOfInputs++){
				while(true){
					int temp = (int) (Math.random() * n);
					if(!inputs.contains(temp)){
						randNum = temp;
						adjMatrix[x][randNum] = true;
						inputs.add(randNum);
						break;
					}
				}
			}
		}

//		// CD4+ T-cell 
//		//adjMatrix[target][source]
//		//0: TBET 
//		//1: IFNG 
//		//2: GATA3 
//		//3: IL2 
//		//4: IL4 
//		//5: RORGT 
//		//6: IL21 
//		//7: FOXP3 
//		//8: TGFB 
//		//9: IL10 
//		//10: BCL6 
//		//11: IL9 
//		//12: IFNGe* 
//		//13: IL2e*
//		//14: IL4e*
//		//15: IL21e*
//		//16: TGFBe*
//		//17: IL10e*
//
//		adjMatrix[0][0] = true;
//		adjMatrix[0][11] = true;
//		adjMatrix[0][4] = true;
//		adjMatrix[0][6] = true;
//		adjMatrix[0][1] = true;
//		adjMatrix[0][2] = true;
//		adjMatrix[0][10] = true;
//
//		adjMatrix[1][8] = true;
//		adjMatrix[1][0] = true;
//		adjMatrix[1][11] = true;
//		adjMatrix[1][4] = true;
//		adjMatrix[1][6] = true;
//		adjMatrix[1][9] = true;
//		adjMatrix[1][12] = true;
//		adjMatrix[1][1] = true;
//		adjMatrix[1][2] = true;
//		adjMatrix[1][10] = true;
//
//		adjMatrix[2][8] = true;
//		adjMatrix[2][0] = true;
//		adjMatrix[2][4] = true;
//		adjMatrix[2][6] = true;
//		adjMatrix[2][3] = true;
//		adjMatrix[2][1] = true;
//		adjMatrix[2][2] = true;
//		adjMatrix[2][10] = true;
//
//		adjMatrix[3][13] = true;
//		adjMatrix[3][6] = true; 
//		adjMatrix[3][3] = true;
//		adjMatrix[3][9] = true; 
//		adjMatrix[3][1] = true;
//		adjMatrix[3][7] = true;
//
//		adjMatrix[4][0] = true;
//		adjMatrix[4][14] = true;
//		adjMatrix[4][4] = true;
//		adjMatrix[4][6] = true;
//		adjMatrix[4][3] = true;
//		adjMatrix[4][1] = true;
//		adjMatrix[4][2] = true;
//
//		adjMatrix[5][8] = true;
//		adjMatrix[5][0] = true;
//		adjMatrix[5][6] = true;
//		adjMatrix[5][2] = true;
//		adjMatrix[5][7] = true;
//		adjMatrix[5][10] = true;
//
//		adjMatrix[6][5] = true;
//		adjMatrix[6][11] = true;
//		adjMatrix[6][4] = true;
//		adjMatrix[6][15] = true;
//		adjMatrix[6][6] = true;
//		adjMatrix[6][3] = true;
//		adjMatrix[6][9] = true;
//		adjMatrix[6][1] = true;
//		adjMatrix[6][10] = true;
//
//		adjMatrix[7][8] = true;
//		adjMatrix[7][5] = true;
//		adjMatrix[7][6] = true;
//		adjMatrix[7][3] = true;
//		adjMatrix[7][7] = true;
//
//		adjMatrix[8][16] = true;
//		adjMatrix[8][8] = true;
//		adjMatrix[8][6] = true;
//		adjMatrix[8][7] = true;
//
//		adjMatrix[9][8] = true;
//		adjMatrix[9][6] = true;
//		adjMatrix[9][17] = true;
//		adjMatrix[9][9] = true;
//		adjMatrix[9][1] = true;
//		adjMatrix[9][2] = true;
//
//		adjMatrix[10][8] = true;
//		adjMatrix[10][0] = true;
//		adjMatrix[10][6] = true;
//		adjMatrix[10][3] = true;
//
//		adjMatrix[11][8] = true;
//		adjMatrix[11][4] = true;
//		adjMatrix[11][6] = true;
//		adjMatrix[11][3] = true;
//		adjMatrix[11][9] = true;
//		adjMatrix[11][1] = true;
//
//		adjMatrix[12][12] = true;
//		adjMatrix[13][13] = true;
//		adjMatrix[14][14] = true;
//		adjMatrix[15][15] = true;
//		adjMatrix[16][16] = true;
//		adjMatrix[17][17] = true;
		
		// state transition table at time T
		if(n <= nodeLimit){
			for (int i = 0; i < nPow; i++) {
				int np = nPow;
				int p = n - 1;
				int a = 0;

				for (int j = 0; j < n; j++) {
					np /= 2;
					nodeVal[i][j] = (i & np) != 0;
					a += (nodeVal[i][j] ? 1 : 0) * (int) Math.pow(2, p);
					p--;
				}

				decNode[i][0] = a;
			}
		}

		// random Boolean functions
		for(int u = 0; u < n; u++){
			for(int kPowIDX = 0; kPowIDX < kPow; kPowIDX++){
				int temp = (int) (Math.random() * 2);
				if(temp == 1){
					booleFunc[kPowIDX][u] = true;
				}
				else{
					booleFunc[kPowIDX][u] = false;
				}
			}
		}
		
//		// CD4+ T-cell 
//		booleFunc[0] = new Vector();
//		for(int idx15 = 0; idx15 < 128; idx15++){
//			if(idx15 == 1 || idx15 == 16 || idx15 == 17){
//				booleFunc[0].addElement(true);
//			}
//			else{
//				booleFunc[0].addElement(false);
//			}
//		}
//
//		booleFunc[1] = new Vector();
//		for(int idx3 = 0; idx3 < 1024; idx3++){
//			if(idx3 == 2 || idx3 == 64 || idx3 == 65 || idx3 == 66 || idx3 == 67 || idx3 == 128 || idx3 == 130
//					|| idx3 == 192 || idx3 == 193 || idx3 == 194 || idx3 == 195	|| idx3 == 320 || idx3 == 321 || idx3 == 322
//					|| idx3 == 323 || idx3 == 448 || idx3 == 449 || idx3 == 450 || idx3 == 451 || idx3 == 576 || idx3 == 577
//					|| idx3 == 578 || idx3 == 579 || idx3 == 704 || idx3 == 705 || idx3 == 706 || idx3 == 707 || idx3 == 832
//					|| idx3 == 833 || idx3 == 834 || idx3 == 835 || idx3 == 960 || idx3 == 961 || idx3 == 962 || idx3 == 963){
//				booleFunc[1].addElement(true);
//			}
//			else{
//				booleFunc[1].addElement(false);
//			}
//		}
//
//		booleFunc[2] = new Vector();
//		for(int idx2 = 0; idx2 < 256; idx2++){
//			if(idx2 == 20 || idx2 == 64 || idx2 == 68 || idx2 == 80 || idx2 == 84){
//				booleFunc[2].addElement(true);
//			}
//			else{
//				booleFunc[2].addElement(false);
//			}
//		}
//
//		booleFunc[3] = new Vector();
//		for(int idx7 = 0; idx7 < 64; idx7++){
//			if(idx7 == 1 || idx7 == 4 || idx7 == 5 || idx7 == 33 || idx7 == 37){
//				booleFunc[3].addElement(true);
//			}
//			else{
//				booleFunc[3].addElement(false);
//			}
//		}
//
//		booleFunc[4] = new Vector();
//		for(int idx11 = 0; idx11 < 128; idx11++){
//			if(idx11 == 2 || idx11 == 3 || idx11 == 6 || idx11 == 7 || idx11 == 18 || idx11 == 19 || idx11 == 22 || idx11 == 23
//					|| idx11 == 66 || idx11 == 67 || idx11 == 68 || idx11 == 70 || idx11 == 71 || idx11 == 80 || idx11 == 82
//					|| idx11 == 83 || idx11 == 84 || idx11 == 86 || idx11 == 87){
//				booleFunc[4].addElement(true);
//			}
//			else{
//				booleFunc[4].addElement(false);
//			}
//		}
//
//		booleFunc[5] = new Vector();
//		for(int idx14 = 0; idx14 < 64; idx14++){
//			if(idx14 == 5){
//				booleFunc[5].addElement(true);
//			}
//			else{
//				booleFunc[5].addElement(false);
//			}
//		}
//
//		booleFunc[6] = new Vector();
//		for(int idx8 = 0; idx8 < 512; idx8++){
//			if(idx8 == 1 || idx8 == 8 || idx8 == 9 || idx8 == 16 || idx8 == 17 || idx8 == 24 || idx8 == 25 || idx8 == 256
//					|| idx8 == 257 || idx8 == 264 || idx8 == 265 || idx8 == 272 || idx8 == 273 || idx8 == 280 || idx8 == 281){
//				booleFunc[6].addElement(true);
//			}
//			else{
//				booleFunc[6].addElement(false);
//			}
//		}
//
//		booleFunc[7] = new Vector();
//		for(int idx1 = 0; idx1 < 32; idx1++){
//			if(idx1 == 9 || idx1 == 24 || idx1 == 25){
//				booleFunc[7].addElement(true);
//			}
//			else{
//				booleFunc[7].addElement(false);
//			}
//		}
//
//		booleFunc[8] = new Vector();
//		for(int idx16 = 0; idx16 < 16; idx16++){
//			if(idx16 == 1 || idx16 == 2 || idx16 == 3 || idx16 == 5 || idx16 == 7 || idx16 == 8 || idx16 == 9 || idx16 == 10
//					|| idx16 == 11 || idx16 == 13 || idx16 == 15){
//				booleFunc[8].addElement(true);
//			}
//			else{
//				booleFunc[8].addElement(false);
//			}
//		}
//
//		booleFunc[9] = new Vector();
//		for(int idx5 = 0; idx5 < 64; idx5++){
//			if(idx5 == 4 || idx5 == 5 || idx5 == 6 || idx5 == 7 || idx5 == 9 || idx5 == 10 || idx5 == 11 || idx5 == 12 || idx5 == 13
//					|| idx5 == 14 || idx5 == 15 || idx5 == 20 || idx5 == 21 || idx5 == 22 || idx5 == 23 || idx5 == 24 || idx5 == 25
//					|| idx5 == 26 || idx5 == 27 || idx5 == 28 || idx5 == 29 || idx5 == 30 || idx5 == 31 || idx5 == 36 || idx5 == 37
//					|| idx5 == 38 || idx5 == 39 || idx5 == 40 || idx5 == 41 || idx5 == 42 || idx5 == 43 || idx5 == 44 || idx5 == 45
//					|| idx5 == 46 || idx5 == 47 || idx5 == 52 || idx5 == 53 || idx5 == 54 || idx5 == 55 || idx5 == 56 || idx5 == 57
//					|| idx5 == 58 || idx5 == 59 || idx5 == 60 || idx5 == 61 || idx5 == 62 || idx5 == 63){
//				booleFunc[9].addElement(true);
//			}
//			else{
//				booleFunc[9].addElement(false);
//			}
//		}
//
//		booleFunc[10] = new Vector();
//		for(int idx0 = 0; idx0 < 16; idx0++){
//			if(idx0 == 4){
//				booleFunc[10].addElement(true);
//			}
//			else{
//				booleFunc[10].addElement(false);
//			}
//		}
//
//		booleFunc[11] = new Vector();
//		for(int idx13 = 0; idx13 < 64; idx13++){
//			if(idx13 == 3 || idx13 == 11 || idx13 == 19 || idx13 == 26 || idx13 == 27){
//				booleFunc[11].addElement(true);
//			}
//			else{
//				booleFunc[11].addElement(false);
//			}
//		}
//
//		booleFunc[12] = new Vector();
//		for(int idx = 0; idx < 2; idx++){
//			if(idx == 0){
//				booleFunc[12].addElement(false);
//			}
//			else{
//				booleFunc[12].addElement(true);
//			}
//		}
//
//		booleFunc[13] = new Vector();
//		for(int idx = 0; idx < 2; idx++){
//			if(idx == 0){
//				booleFunc[13].addElement(false);
//			}
//			else{
//				booleFunc[13].addElement(true);
//			}
//		}
//
//		booleFunc[14] = new Vector();
//		for(int idx = 0; idx < 2; idx++){
//			if(idx == 0){
//				booleFunc[14].addElement(false);
//			}
//			else{
//				booleFunc[14].addElement(true);
//			}
//		}
//
//		booleFunc[15] = new Vector();
//		for(int idx = 0; idx < 2; idx++){
//			if(idx == 0){
//				booleFunc[15].addElement(false);
//			}
//			else{
//				booleFunc[15].addElement(true);
//			}
//		}
//
//		booleFunc[16] = new Vector();
//		for(int idx = 0; idx < 2; idx++){
//			if(idx == 0){
//				booleFunc[16].addElement(false);
//			}
//			else{
//				booleFunc[16].addElement(true);
//			}
//		}
//
//		booleFunc[17] = new Vector();
//		for(int idx = 0; idx < 2; idx++){
//			if(idx == 0){
//				booleFunc[17].addElement(false);
//			}
//			else{
//				booleFunc[17].addElement(true);
//			}
//		}	
	}

	public void makeBooleanNetwork() 
	{

		for(int i = 0; i < n; i++){
			appliedNode[i] = new Vector();
			for(int j = 0; j < n; j++){
				if(adjMatrix[i][j] == true){
					appliedNode[i].addElement(j);
				}
			}
		}
		
//		// CD4+ T-cell
//		for(int i = 0; i < n; i++){
//			appliedNode[i] = new Vector();
//		}
//		appliedNode[0].addElement(10);
//		appliedNode[0].addElement(2);
//		appliedNode[0].addElement(1);
//		appliedNode[0].addElement(6);
//		appliedNode[0].addElement(4);
//		appliedNode[0].addElement(11);
//		appliedNode[0].addElement(0);
//
//		appliedNode[1].addElement(10);
//		appliedNode[1].addElement(2);
//		appliedNode[1].addElement(1);
//		appliedNode[1].addElement(12);
//		appliedNode[1].addElement(9);
//		appliedNode[1].addElement(6);
//		appliedNode[1].addElement(4);
//		appliedNode[1].addElement(11);
//		appliedNode[1].addElement(0);
//		appliedNode[1].addElement(8);
//
//		appliedNode[2].addElement(10);
//		appliedNode[2].addElement(2);
//		appliedNode[2].addElement(1);
//		appliedNode[2].addElement(3);
//		appliedNode[2].addElement(6);
//		appliedNode[2].addElement(4);
//		appliedNode[2].addElement(0);
//		appliedNode[2].addElement(8);
//
//		appliedNode[3].addElement(7);
//		appliedNode[3].addElement(1);
//		appliedNode[3].addElement(9);
//		appliedNode[3].addElement(3);
//		appliedNode[3].addElement(6);
//		appliedNode[3].addElement(13);
//
//		appliedNode[4].addElement(2);
//		appliedNode[4].addElement(1);
//		appliedNode[4].addElement(3);
//		appliedNode[4].addElement(6);
//		appliedNode[4].addElement(4);
//		appliedNode[4].addElement(14);
//		appliedNode[4].addElement(0);
//
//		appliedNode[5].addElement(10);
//		appliedNode[5].addElement(7);
//		appliedNode[5].addElement(2);
//		appliedNode[5].addElement(6);
//		appliedNode[5].addElement(0);
//		appliedNode[5].addElement(8);
//
//		appliedNode[6].addElement(10);
//		appliedNode[6].addElement(1);
//		appliedNode[6].addElement(9);
//		appliedNode[6].addElement(3);
//		appliedNode[6].addElement(6);
//		appliedNode[6].addElement(15);
//		appliedNode[6].addElement(4);
//		appliedNode[6].addElement(11);
//		appliedNode[6].addElement(5);
//
//		appliedNode[7].addElement(7);
//		appliedNode[7].addElement(3);
//		appliedNode[7].addElement(6);
//		appliedNode[7].addElement(5);
//		appliedNode[7].addElement(8);
//
//		appliedNode[8].addElement(7);
//		appliedNode[8].addElement(6);
//		appliedNode[8].addElement(8);
//		appliedNode[8].addElement(16);
//
//		appliedNode[9].addElement(2);
//		appliedNode[9].addElement(1);
//		appliedNode[9].addElement(9);
//		appliedNode[9].addElement(17);
//		appliedNode[9].addElement(6);
//		appliedNode[9].addElement(8);
//
//		appliedNode[10].addElement(3);
//		appliedNode[10].addElement(6);
//		appliedNode[10].addElement(0);
//		appliedNode[10].addElement(8);
//
//		appliedNode[11].addElement(1);
//		appliedNode[11].addElement(9);
//		appliedNode[11].addElement(3);
//		appliedNode[11].addElement(6);
//		appliedNode[11].addElement(4);
//		appliedNode[11].addElement(8);
//
//		appliedNode[12].addElement(12);
//		appliedNode[13].addElement(13);
//		appliedNode[14].addElement(14);
//		appliedNode[15].addElement(15);
//		appliedNode[16].addElement(16);
//		appliedNode[17].addElement(17);
		
		if(n <= nodeLimit){
			// state transition table at time T+1
			newNodeVal = new boolean[nPow][n];
			newDecNode = new int[nPow][1];
			for(int i = 0; i < n; i++){
				for(int j = 0; j < nPow; j++){

					// random Boolean networks
					int a = 0;
					for(int p = 0; p < appliedNode[i].size(); p++){
						int temp = appliedNode[i].get(appliedNode[i].size() - p - 1);
						a += (nodeVal[j][temp]?1:0) * (int) Math.pow(2, p);
					}
					boolean value = booleFunc[a][i];
//					boolean value = booleFunc[i].get(a); // CD4+ T-cell

					newNodeVal[j][i] = value;
				}
			}

			for(int i = 0; i < nPow; i++){
				int a = 0;
				int p = n-1;
				for(int j = 0; j < n; j++){
					a += (newNodeVal[i][j]?1:0) * (int) Math.pow(2,p);
					p--;
				}
				newDecNode[i][0] = a;
			}
		}	
	}

	public void calComplexity(){


		// 2D array to store random initial states 
		int[][] randomInitial = new int[powLimit][1];

		// initialization for randomInitial array
		for(int ii = 0; ii < powLimit; ii++){
			randomInitial[ii] = null;
		}	

		// choosing random initial states
		for(int i = 0; i < powLimit; i++){

			int[] tempRandom = new int[n];

			for(int j = 0; j < n; j++){
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

		// // tracking state transitions until T
		int[][] configurationTrack = new int[(simulTime_T*2)+1][1]; 
		int[][] configurationTrack2 = new int[(simulTime_T*2)+1][1];

		// before adding perturbations
		int[] tempConfigBinary = new int[n];	
		int[] tempNewConfigBinary = new int[n];

		// after adding perturbations
		int[] tempConfigBinary2 = new int[n];	
		int[] tempNewConfigBinary2 = new int[n];

		for(int conIdx = 0; conIdx < powLimit; conIdx++){

			// initialization of configurationTrack array
			for(int ii = 0; ii < simulTime_T*2; ii++){
				configurationTrack[ii] = null;
				configurationTrack2[ii] = null;	
			}

			configurationTrack[0] = randomInitial[conIdx];
			configurationTrack2[0] = randomInitial[conIdx];	// same random initial states

			for(int init1 = 0; init1 < n; init1++){	// initialization for tempConfigBinary array
				tempConfigBinary2[init1] = randomInitial[conIdx][init1];
			}			

			for(int init2 = 0; init2 < n; init2++){	// initialization for tempNewConfigBinary array
				tempNewConfigBinary[init2] = -1;
				tempNewConfigBinary2[init2] = -1;
			}

			int time = 0;

			while(time < simulTime_T*2){

				if(n <= nodeLimit){
					int tempDecVal = 0;
					int pp = n - 1;
					for(int bitIdx = 0; bitIdx < n; bitIdx++){
						tempDecVal = tempDecVal + tempConfigBinary[bitIdx] * (int) Math.pow(2,pp);
						pp--;
					}

					int tempDecVal2 = 0;
					int pp2 = n - 1;
					for(int bitIdx2 = 0; bitIdx2 < n; bitIdx2++){
						tempDecVal2 = tempDecVal2 + tempConfigBinary2[bitIdx2] * (int) Math.pow(2,pp2);
						pp2--;
					}

					for(int idx = 0; idx < n; idx++){
						tempNewConfigBinary[idx] = (newNodeVal[tempDecVal][idx]?1:0);
						tempNewConfigBinary2[idx] = (newNodeVal[tempDecVal2][idx]?1:0);
					}
				}
				else{ 
					for(int i = 0; i < n; i++){
						// random Boolean networks
						int a = 0;
						int a2 = 0;
						int len = appliedNode[i].size()-1;
						for(int p = 0; p < appliedNode[i].size() ; p++){
							int idx = appliedNode[i].get(len);
							a += (tempConfigBinary[idx]) * (int) Math.pow(2, p);
							a2 += (tempConfigBinary2[idx]) * (int) Math.pow(2, p);
							len--;
						}

						boolean value = booleFunc[a][i];
//						boolean value = booleFunc[i].get(a); // CD4+ T-cell
						
						if(value == true){
							tempNewConfigBinary[i] = 1;
						}
						else{
							tempNewConfigBinary[i] = 0;
						}

						boolean value2 = booleFunc[a2][i]; 
//						boolean value2 = booleFunc[i].get(a2); // CD4+ T-cell 
						
						if(value2 == true){
							tempNewConfigBinary2[i] = 1;
						}
						else{
							tempNewConfigBinary2[i] = 0;
						}
					}
				}

				int[] temp = new int[n];
				int[] temp2 = new int[n];

				int1DArrayCopier(temp,tempNewConfigBinary, n);
				int1DArrayCopier(temp2,tempNewConfigBinary2, n);
				time++;

				// CD4+ T-cell: the states of the nodes without inputs continue to be maintained as the values of initial states 
//				temp2[12] = configurationTrack2[0][12];
//				temp2[13] = configurationTrack2[0][13];
//				temp2[14] = configurationTrack2[0][14];
//				temp2[15] = configurationTrack2[0][15];
//				temp2[16] = configurationTrack2[0][16];
//				temp2[17] = configurationTrack2[0][17];
				
				configurationTrack[time] = temp;	
				if(time%perPeriod_O == 0){	
					// choosing perturbed nodes
					Vector<Integer> randNodeIDX4Per = new Vector<Integer>();
					for(int p = 0; p < perSize_X; p++){ 

						int tempRandIDX = (int)(Math.random()*n);
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
				int1DArrayCopier(tempConfigBinary,temp, n);	 
				int1DArrayCopier(tempConfigBinary2,temp2, n);
			}

			// measuring complexity based on Shannon's entropy 
			double totalE = 0;
			double totalE2 = 0;
			double initialComplexity = 0;
			double finalComplexity = 0;
			for(int nIDX = 0; nIDX < n; nIDX++){
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
			initialComplexity = 4*((double)totalE/n)*(1-((double)totalE/n));
			finalComplexity = 4*((double)totalE2/n)*(1-((double)totalE2/n));
			initialComplexityStorage.addElement(initialComplexity);
			finalComplexityStorage.addElement(finalComplexity);
			delSigma.addElement(finalComplexity-initialComplexity);
		}
	}

	public void calFragility(){

		double degreeOfPer = (perSize_X*((double) simulTime_T/perPeriod_O))/((double) n*simulTime_T);
		double tempFragility = -2;
		for(int cv = 0; cv < delSigma.size(); cv++){
			tempFragility = -delSigma.get(cv) * degreeOfPer;
			fragility.addElement(tempFragility);
		}
	}

	public void individualNetwork(){
		makeStructure();
		makeBooleanNetwork();
	}

	public void MeasureAntifragile(){
		calComplexity(); 
		calFragility();
	}

	double log(double base, double num)
	{
		return Math.log(num) / Math.log(base);
	}

	public void intVectorCopier(Vector<Integer> targetVector, Vector<Integer> sourceVector)
	{
		Integer x;
		for(int i = 0; i < sourceVector.size(); i++)
		{
			x = new Integer(sourceVector.get(i).intValue());
			targetVector.add(x);
		}
	}

	public void intVectorArrayCopier(Vector<Integer>[] targetVectorArray, Vector<Integer>[] sourceVectorArray)
	{
		for(int i = 0; i < n; i++){
			targetVectorArray[i] = new Vector<Integer>();
			for(int j = 0; j < sourceVectorArray[i].size(); j++){
				targetVectorArray[i].add(sourceVectorArray[i].get(j));
			}

		}
	}

	public void boolean2DArrayCopier(boolean targetArray[][], boolean sourceArray[][], int x, int y)
	{
		for(int i = 0; i < x; i++)
		{
			for(int j = 0; j < y; j++)
			{
				targetArray[i][j] = sourceArray[i][j];
			}
		}
	}

	public void int2DArrayCopier(int targetArray[][], int sourceArray[][], int x,int y)
	{
		for(int i = 0; i < x; i++)
		{
			for(int j = 0; j < y; j++)
			{
				targetArray[i][j] = sourceArray[i][j];
			}
		}
	}

	public void int1DArrayCopier(int targetArray[], int sourceArray[], int x)
	{
		for(int i = 0; i < x; i++)
		{
			targetArray[i] = sourceArray[i];
		}
	}

	public void double1DArrayCopier(double targetArray[], double sourceArray[], int x)
	{
		for(int i = 0; i < x; i++)
		{
			targetArray[i] = sourceArray[i];
		}
	}

}