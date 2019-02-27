import java.io.*;
import java.util.Vector;

public class Driver {

	static long start = System.currentTimeMillis();
	BufferedWriter output;

	public Driver() throws IOException{
		output = new BufferedWriter(new FileWriter("./output.txt"));
	}

	public Driver(String outFilename) throws IOException{
		output = new BufferedWriter(new FileWriter(outFilename));
	}

	public void process() throws IOException{

		//single-layer RBNs
//		int iteration = 100;
//		int infoFlag = 0;
//		Vector<Double> finalComplexityVector = new Vector<Double>();
//		Vector<Double> initialComplexityVector = new Vector<Double>();
//		Vector<Double> fragilityVector = new Vector<Double>();
//		int numOfAntiFragile = 0;
//		double sumOfFragility = 0;
//
//		// perturbation size X		
//		Vector<Integer> perSizeVector = new Vector<Integer>();
//		for(int psIDX = 1; psIDX < 163; psIDX++){
//			perSizeVector.add(psIDX);
//		}
//		
//		// perturbation frequency O
//		Vector<Integer> perPeriodVector = new Vector<Integer>();		
//		perPeriodVector.add(1);
//		perPeriodVector.add(2);
//		perPeriodVector.add(3);
//		perPeriodVector.add(4);
//		perPeriodVector.add(5);
//
//		for(int perPeriodIDX = 0; perPeriodIDX < perPeriodVector.size(); perPeriodIDX++){
//			for(int perSizeIDX = 0; perSizeIDX < perSizeVector.size(); perSizeIDX++){
//				
//				infoFlag = 0;
//				finalComplexityVector = new Vector<Double>();
//				initialComplexityVector = new Vector<Double>();
//				fragilityVector = new Vector<Double>();
//				numOfAntiFragile = 0;
//				sumOfFragility = 0;
//				
//				for(int i = 0; i < iteration; i++){
//					
//					System.out.println("==iteration = "+i);
//					NKBoolean nk = new NKBoolean();
//					nk.individualNetwork();
//					nk.perPeriod_O = perPeriodVector.get(perPeriodIDX);
//					nk.perSize_X = perSizeVector.get(perSizeIDX);
//					nk.MeasureAntifragile();
//					
//					if(infoFlag == 0){
//						String s1 = "n = "+nk.n +", k = "+nk.k+", T = "+nk.simulTime_T+", X = "+nk.perSize_X+", O = "+nk.perPeriod_O+", numberOfInitials = "+nk.powLimit+", iteration = "+iteration;
//						output.write(s1); output.newLine();
//						infoFlag = 1;
//					}
//					String s2 = "== iteration = "+i;
//					output.write(s2); output.newLine();
//
//					String s3 = "finalComplexity, ";
//					output.write(s3); 
//					double temp1 = 0;
//					for(int p1 = 0; p1 < nk.finalComplexityStorage.size(); p1++){
//						temp1 = temp1 + nk.finalComplexityStorage.get(p1);
//						String s4 = nk.finalComplexityStorage.get(p1)+", ";
//						output.write(s4); 
//					}
//					double avgFinalComplexity = (double) temp1/nk.finalComplexityStorage.size();
//					finalComplexityVector.add(avgFinalComplexity);
//					output.newLine();
//
//					String s5 = "initialComplexity, ";
//					output.write(s5);
//					double temp2 = 0;
//					for(int p2 = 0; p2 < nk.initialComplexityStorage.size(); p2++){
//						temp2 = temp2 + nk.initialComplexityStorage.get(p2);
//						String s6 = nk.initialComplexityStorage.get(p2)+", ";
//						output.write(s6); 
//					}
//					double avgInitialComplexity = (double) temp2/nk.initialComplexityStorage.size();
//					initialComplexityVector.add(avgInitialComplexity);
//					output.newLine();
//
//					String s7 = "fragility, ";
//					output.write(s7); 
//					double temp3 = 0;
//					for(int p3 = 0; p3 < nk.fragility.size(); p3++){
//						temp3 = temp3 + nk.fragility.get(p3);
//						String s8 = nk.fragility.get(p3)+", ";
//						output.write(s8); 
//					}
//					double avgFragility = (double) temp3/nk.fragility.size();
//					fragilityVector.add(avgFragility);
//					output.newLine();
//
//					if(avgFragility < 0){
//						numOfAntiFragile++;
//					}
//					sumOfFragility = sumOfFragility + avgFragility;
//				}
//
//				double meanFragility = (double) sumOfFragility/iteration;
//				double standardDeviation = 0;
//				double stdFragility = 0;
//				for(int idx = 0; idx < iteration; idx++) {
//					standardDeviation += Math.pow(fragilityVector.get(idx) - meanFragility, 2);
//				}
//				stdFragility = Math.sqrt((double) standardDeviation/(iteration-1));
//
//				String s8 = "== summary ";
//				output.write(s8); output.newLine();
//
//				String s9 = "avgFinalComplexity = "+finalComplexityVector;
//				output.write(s9); output.newLine();
//
//				String s10 = "avgInitialComplexity = "+initialComplexityVector;
//				output.write(s10); output.newLine();
//
//				String s11 = "avgFragility = "+fragilityVector;
//				output.write(s11); output.newLine();
//
//				String s12 = "numOfAntifragile = "+numOfAntiFragile;
//				output.write(s12); output.newLine();
//
//				String s13 = "avgAvgFragility = " + meanFragility;
//				output.write(s13); output.newLine();
//
//				String s14 = "stdAvgFragility = " + stdFragility;
//				output.write(s14); output.newLine(); output.newLine();
//			}
//		}


		
		// multilayer RBNs
		int iteration = 100;
		int infoFlag = 0;
		Vector<Double> finalComplexityVector = new Vector<Double>();
		Vector<Double> initialComplexityVector = new Vector<Double>();
		Vector<Double> fragilityVector = new Vector<Double>();
		int numOfAntiFragile = 0;
		double sumOfFragility = 0;
		
		// perturbation size X
		Vector<Integer> perSizeVector = new Vector<Integer>();
		for(int psIDX = 1; psIDX < 163 ; psIDX++){
			perSizeVector.add(psIDX);
		}
		
		// perturbation frequency O
		Vector<Integer> perPeriodVector = new Vector<Integer>();
		perPeriodVector.add(1);
		perPeriodVector.add(2);
		perPeriodVector.add(3);
		perPeriodVector.add(4);
		perPeriodVector.add(5);

		for(int perPeriodIDX = 0; perPeriodIDX < perPeriodVector.size(); perPeriodIDX++){
			for(int perSizeIDX = 0; perSizeIDX < perSizeVector.size(); perSizeIDX++){
				
				infoFlag = 0;
				finalComplexityVector = new Vector<Double>();
				initialComplexityVector = new Vector<Double>();
				fragilityVector = new Vector<Double>();
				numOfAntiFragile = 0;
				sumOfFragility = 0;
			
				for(int i = 0; i < iteration; i++){
					System.out.println("==iteration = "+i);
					HierarchicalNKBoolean hNK = new HierarchicalNKBoolean();					
					hNK.perPeriod_O = perPeriodVector.get(perPeriodIDX);
					hNK.perSize_X = perSizeVector.get(perSizeIDX);
					hNK.hierarchicalNetwork();
					
					if(infoFlag == 0){
						String s1 = "n = "+hNK.n +", k = "+hNK.k+", interNode = "+hNK.interN+", interLink = "+hNK.interLink+", T = "+hNK.simulTime_T+
								", X = "+hNK.perSize_X+", O = "+hNK.perPeriod_O+", numberOfInitials = "+hNK.hPowLimit+
								", numberOfCommunicatingNodes = "+hNK.communicatingNodeIDXVector.size()+", iteration = "+iteration;
						output.write(s1); output.newLine();
						infoFlag = 1;
					}
					String s2 = "== iteration = "+i;
					output.write(s2); output.newLine();

					String s3 = "finalComplexity, ";
					output.write(s3); 
					double temp1 = 0;
					for(int p1 = 0; p1 < hNK.finalComplexityStorage.size(); p1++){
						temp1 = temp1 + hNK.finalComplexityStorage.get(p1);
						String s4 = hNK.finalComplexityStorage.get(p1)+", ";
						output.write(s4); 
					}
					double avgFinalComplexity = (double) temp1/hNK.finalComplexityStorage.size();
					finalComplexityVector.add(avgFinalComplexity);
					output.newLine();

					String s5 = "initialComplexity, ";
					output.write(s5);
					double temp2 = 0;
					for(int p2 = 0; p2 < hNK.initialComplexityStorage.size(); p2++){
						temp2 = temp2 + hNK.initialComplexityStorage.get(p2);
						String s6 = hNK.initialComplexityStorage.get(p2)+", ";
						output.write(s6); 
					}
					double avgInitialComplexity = (double) temp2/hNK.initialComplexityStorage.size();
					initialComplexityVector.add(avgInitialComplexity);
					output.newLine();

					String s7 = "fragility, ";
					output.write(s7); 
					double temp3 = 0;
					for(int p3 = 0; p3 < hNK.fragility.size(); p3++){
						temp3 = temp3 + hNK.fragility.get(p3);
						String s8 = hNK.fragility.get(p3)+", ";
						output.write(s8); 
					}
					double avgFragility = (double) temp3/hNK.fragility.size();
					fragilityVector.add(avgFragility);
					output.newLine();

					if(avgFragility < 0){
						numOfAntiFragile++;
					}

					sumOfFragility = sumOfFragility + avgFragility;
				}
				
				double meanFragility = (double) sumOfFragility/iteration;
				double standardDeviation = 0;
				double stdFragility = 0;
				for(int idx = 0; idx < iteration; idx++) {
					standardDeviation += Math.pow(fragilityVector.get(idx) - meanFragility, 2);
				}
				stdFragility = Math.sqrt((double) standardDeviation/(iteration-1));

				String s8 = "== summary ";
				output.write(s8); output.newLine();

				String s9 = "avgFinalComplexity = "+finalComplexityVector;
				output.write(s9); output.newLine();

				String s10 = "avgInitialComplexity = "+initialComplexityVector;
				output.write(s10); output.newLine();

				String s11 = "avgFragility = "+fragilityVector;
				output.write(s11); output.newLine();

				String s12 = "numOfAntifragile = "+numOfAntiFragile;
				output.write(s12); output.newLine();

				String s13 = "avgAvgFragility = " + (double) sumOfFragility/iteration;
				output.write(s13); output.newLine();

				String s14 = "stdAvgFragility = " + stdFragility;
				output.write(s14); output.newLine(); output.newLine();
			}
		}
		
		

		long end = System.currentTimeMillis();		

		String s2 = end-start+" ms passed"; //minutes: *(1/60000)
		output.write(s2); output.newLine();
		output.close();
	}

	public static void main(String[] args) throws IOException{

		for(int j = 0; j < 1; j++){

			File dir = new File("./");
			File[] files = dir.listFiles();

			String newFilename = null;
			if(files == null){
				throw new IOException("aaa");
			}
			else{

				int maxIdx = -1;
				for( int w = 0; w < files.length; w++){
					String exFile = files[w].getName();

					if(exFile.endsWith(".jar")){
						continue;
					}

					if(!exFile.startsWith("output_")){
						continue;
					}

					if(!exFile.endsWith(".txt")){
						continue;
					}

					int idx1 = "output_".length();
					int idx2 = exFile.lastIndexOf(".txt");
					int currIdx = Integer.parseInt(exFile.substring(idx1, idx2));

					if(currIdx > maxIdx){
						maxIdx = currIdx;
					}
				}

				if(maxIdx == -1){
					newFilename = "output_1.txt";
				}
				else{
					newFilename = "output_"+ String.valueOf(++maxIdx) + ".txt"; 
				}
			}

			Driver mainDriver = new Driver(newFilename);
			mainDriver.process();
		}
		System.out.println("--END--");
	}
}
