package Model;

import java.util.*;
import Data.*;

public class PositionDistributionStatistic extends LogLikelihoodRatio {
	public int bestWindowSize;
	public int bestOffset;
	public double caseUnormLogLikelihood;
	public double controlUnormLogLikelihood;
	public double caseNormConstant;
	public double controlNormConstant;
	public double allUnormLogLikelihood;
	
	public PositionDistributionStatistic(){
		this.bestOffset = -1;
		this.bestWindowSize = -1;
		this.allUnormLogLikelihood = 0.0;
		this.caseNormConstant = 0.0;
		this.controlUnormLogLikelihood = 0.0;
		this.caseUnormLogLikelihood = 0.0;
		this.controlNormConstant = 0.0;
		this.LLAlternative = 1.0;
		this.LLNull = 1.0;
	}
	
	//compute position distribution statistic for one genomic unit
	public PositionDistributionStatistic(LinkedList<Double> affectedCounts, LinkedList<Double> unaffectedCounts, LinkedList<Variant> variants, GenomicUnit unit, Setting setting){
		this.LLAlternative = 0;
		this.LLNull = 0;
		this.bestWindowSize = setting.windowLength;
		this.bestOffset = 0;
		
		//the number of subunits defined in the genomic unit
		int numOfSubUnit = 0;
		if(unit.subUnits != null) numOfSubUnit = unit.subUnits.length;
		
		//get pseudo count 
		double zero = pseudoCount(setting.positionScoring);
		
		//***************debug****************
//		if(unit.name.equals("MAP3K6")){
//			System.out.println(unit.name+","+String.valueOf(unit.length)+", #subunit: "+String.valueOf(numOfSubUnit));
//			System.out.println("pseudoCount: "+String.valueOf(zero));
//			
//			Iterator<Variant> variantItr = variants.iterator();
//			Iterator<Double> affectedCountItr = affectedCounts.iterator();
//			Iterator<Double> unaffectedCountItr = unaffectedCounts.iterator();
//			
//			while(variantItr.hasNext()){
//				Variant var = variantItr.next();
//				double affectedCount = affectedCountItr.next();
//				double unaffectedCount = unaffectedCountItr.next();
//				
//				System.out.println("("+String.valueOf(var.pos)+","+String.valueOf(var.preCompPositionScore)+","+
//				String.valueOf(affectedCount)+","+String.valueOf(unaffectedCount)+")");
//			}
//		}
		//***************debug****************
		
		//keep the largest log likelihood ratios
		//one is for all variants (Acc) and the other is for unique variants only
		//the best segmentation is picked according to one of these ratios
		double currBestAccLogRatio = 0.0;
		double currBestLogRatio = 0.0;
		
		//currWindowSize controls the window size used to compute local likelihood.
		//It grows with a factor of 2 in each iteration.
		int currWindowSize = setting.windowLength;
		
		//offsetStep is how many bases the window segmentation shifts every iteration
		int offsetStep = setting.windowOffset;
		
		for(int i = 0;i < setting.windowGrowth; i++){
			
			//the window segmentation is shifted by offsetStep in every iteration
			for(int currOffset = 0;currOffset < currWindowSize;currOffset += offsetStep){
				
				//used in calculation of log likelihood ratio for unique variants only
				double currLogLikelihoodRatio = 0.0;
				double currCaseLogLikelihood = 0.0;
				double currControlLogLikelihood = 0.0;
				double currAllLogLikelihood = 0.0;

				//used in calculation of log likelihood ratio for all variants
				double currAccLogLikelihoodRatio = 0.0;
				double currCaseAccLogLikelihood = 0.0;
				double currControlAccLogLikelihood = 0.0;
				double currAllAccLogLikelihood = 0.0;
				
				//calculate how many windows segment the genomic unit
				int windowLength = (unit.length-currOffset)/currWindowSize+2;
				
				//the genomic unit is segmented by fixed length windows and subunits
				//for unique variants only
				double[] caseWindowScore = new double[windowLength+numOfSubUnit];
				double[] controlWindowScore = new double[windowLength+numOfSubUnit];
				double[] allWindowScore = new double[windowLength+numOfSubUnit];
				
				//for all variants
				double[] caseAccWindowScore = new double[windowLength+numOfSubUnit];
				double[] controlAccWindowScore = new double[windowLength+numOfSubUnit];
				double[] allAccWindowScore = new double[windowLength+numOfSubUnit];
				
				//pseudo count for each window
				for(int j=0;j<(windowLength+numOfSubUnit);j++){
					caseWindowScore[j] = zero;
					controlWindowScore[j] = zero;
					allWindowScore[j] = 2*zero;
					
					caseAccWindowScore[j] = zero;
					controlAccWindowScore[j] = zero;
					allAccWindowScore[j] = 2*zero;
				}
				
				//add scores into proper windows one variant by one variant
				Iterator<Variant> variantItr = variants.iterator();
				Iterator<Double> affectedCountItr = affectedCounts.iterator();
				Iterator<Double> unaffectedCountItr = unaffectedCounts.iterator();
				
				while(variantItr.hasNext()){
					Variant entryVar = variantItr.next();
					double affectedScore = affectedCountItr.next();
					double unaffectedScore = unaffectedCountItr.next();
					
					//window index
					int pos = 0;
					
					if(entryVar.subGenomicUnitId == -1){
						//variant in main region of the genomic unit so calculate the proper window based on its position
						pos = (entryVar.pos-currOffset-1)/currWindowSize;
						if(currOffset != 0) pos++;
					} else {
						//variant in a subunit of the genomic unit so put it in the corresponding window
						pos = windowLength+entryVar.subGenomicUnitId;
					}
					
					//add scores accounting for all variants at this position in the study
					caseAccWindowScore[pos] += affectedScore;
					controlAccWindowScore[pos] += unaffectedScore;
					allAccWindowScore[pos] += (affectedScore+unaffectedScore);
					
					//for the calculation of unique variant score, the outcome of one mutational position
					//is either "present" or "not present"
					double affectedIndicator = (affectedScore > 0 ? entryVar.preCompPositionScore : 0);
					double unaffectedIndicator = (unaffectedScore > 0 ? entryVar.preCompPositionScore : 0);
					
					//add scores accounting for unique variants only at this position in the study
					caseWindowScore[pos] += affectedIndicator;
					controlWindowScore[pos] += unaffectedIndicator;
					allWindowScore[pos] += (affectedIndicator+unaffectedIndicator);
				}
				
				double caseScoreSum = 0.0;
				double controlScoreSum = 0.0;
				
				double caseAccScoreSum = 0.0;
				double controlAccScoreSum = 0.0;
				
				//compute numerators of likelihoods and the normalization factors
				for(int j=0;j<(windowLength+numOfSubUnit);j++){
					if(caseWindowScore[j] > 0){
						currCaseLogLikelihood += caseWindowScore[j]*Math.log(caseWindowScore[j]);
						caseScoreSum += caseWindowScore[j];
						
						currCaseAccLogLikelihood += caseAccWindowScore[j]*Math.log(caseAccWindowScore[j]);
						caseAccScoreSum += caseAccWindowScore[j];
					}
					if(controlWindowScore[j] > 0){
						currControlLogLikelihood += controlWindowScore[j]*Math.log(controlWindowScore[j]);
						controlScoreSum += controlWindowScore[j];
						
						currControlAccLogLikelihood += controlAccWindowScore[j]*Math.log(controlAccWindowScore[j]);
						controlAccScoreSum += controlAccWindowScore[j];
					}
					if(caseWindowScore[j] > 0 || controlWindowScore[j] > 0){
						currAllLogLikelihood += (caseWindowScore[j]+controlWindowScore[j])*Math.log(caseWindowScore[j]+controlWindowScore[j]);
						currAllAccLogLikelihood += (caseAccWindowScore[j]+controlAccWindowScore[j])*Math.log(caseAccWindowScore[j]+controlAccWindowScore[j]);
					}
				}
				
				//normalized by normalization factors
				if(caseScoreSum != 0) currCaseLogLikelihood -= caseScoreSum*Math.log(caseScoreSum);
				if(controlScoreSum != 0) currControlLogLikelihood -= controlScoreSum*Math.log(controlScoreSum);
				currAllLogLikelihood -= (caseScoreSum+controlScoreSum)*Math.log(caseScoreSum+controlScoreSum);
				
				currLogLikelihoodRatio =currCaseLogLikelihood+currControlLogLikelihood-currAllLogLikelihood;
				
				double currCaseUnormAccLogLikelihood = currCaseAccLogLikelihood;
				double currControlUnormAccLogLikelihood = currControlAccLogLikelihood;
				double currAllUnormAccLogLikelihood = currAllAccLogLikelihood;
				
				if(caseAccScoreSum != 0) currCaseAccLogLikelihood -= caseAccScoreSum*Math.log(caseAccScoreSum);
				if(controlAccScoreSum != 0) currControlAccLogLikelihood -= controlAccScoreSum*Math.log(controlAccScoreSum);
				currAllAccLogLikelihood -= (caseAccScoreSum+controlAccScoreSum)*Math.log(caseAccScoreSum+controlAccScoreSum);
				
				currAccLogLikelihoodRatio =currCaseAccLogLikelihood+currControlAccLogLikelihood-currAllAccLogLikelihood;
				
				//see if we get a better window segmentation which generates a greater log likelihood ratio
				boolean needUpdate = false;
				if(setting.bestSegmentationByUniqueVar){
					if(currLogLikelihoodRatio > currBestLogRatio) needUpdate = true;
				} else {
					if(currAccLogLikelihoodRatio > currBestAccLogRatio) needUpdate = true;
				}
				
				if(needUpdate){
					this.LLAlternative = currCaseAccLogLikelihood+currControlAccLogLikelihood;
					this.LLNull = currAllAccLogLikelihood;
					this.bestWindowSize = currWindowSize;
					this.bestOffset = currOffset;
					this.caseUnormLogLikelihood = currCaseUnormAccLogLikelihood;
					this.controlUnormLogLikelihood = currControlUnormAccLogLikelihood;
					this.caseNormConstant = caseAccScoreSum;
					this.controlNormConstant = controlAccScoreSum;
					this.allUnormLogLikelihood = currAllUnormAccLogLikelihood;
					currBestLogRatio = currLogLikelihoodRatio;
					currBestAccLogRatio = currAccLogLikelihoodRatio;
				}
			}
			
			//grows the window by the factor of 2
			currWindowSize *= 2;
			offsetStep *= 2;
		}
	}
	
	double pseudoCount(int scoringMode){
		return 1.0;
	}

}
