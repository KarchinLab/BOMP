package Model;

import java.util.*;
import Data.*;

//This object picks a specific burden score for a genomic unit
//every individual who has a score higher than or equal to the cutoff is considered having a functional damage
//the cutoff is picked by maximizing the log likelihood ratio

public class BurdenHTStatistic extends LogLikelihoodRatio {
	public double aboveThresholdA;
	public double aboveThresholdU;
	public double cutoff;
	
	public BurdenHTStatistic(Map<Integer,Double> sampleId2score, CaseControlStudy study){
		SortedSet<Map.Entry<Integer, Double>> sortedSet = Utility.entriesSortedByValues(sampleId2score);
		
		this.aboveThresholdA = -1;
		this.aboveThresholdU = -1;
		this.LLAlternative = 0;
		this.LLNull = 0;
		
		Iterator<Map.Entry<Integer, Double>> entryItr = sortedSet.iterator();
		Map.Entry<Integer, Double> one = entryItr.next();
		int sampleId = one.getKey();
		double score = one.getValue();
		double currCutOff = score;
		int currAffectedCase = 0;
		int currAffectedControl = 0;
		if(study.sampleId2affected.get(sampleId)) currAffectedCase++;
		else currAffectedControl++;
		
		//================ debug =================
		//Sample testOne = study.sampleId2sample.get(sampleId);
		//System.out.println(testOne.name+": "+String.valueOf(score));
		//================ debug =================
		
		//processing the following if they exist
		while(entryItr.hasNext()){
			one = entryItr.next();
			sampleId = one.getKey();
			score = one.getValue();
			
			//================ debug =================
			//testOne = study.sampleId2sample.get(sampleId);
			//System.out.println(testOne.name+": "+String.valueOf(score));
			//================ debug =================
			
			if(currCutOff == score){
				
				//the next one has the same burden as the current cutoff
				//samples are affected if their burdens are greater than or equal to the cutoff
				
				if(study.sampleId2affected.get(sampleId)) currAffectedCase++;
				else currAffectedControl++;
				
			} else if(currCutOff > score){
				
				//find the boundary of whether the sample is affected or not for the current cutoff
				//compute likelihood ratio and keep it if it is greater than the previous one
				
				//adding pseudo counts
				//reflecting the prior belief of whether this gene is related to the phenotype or not

				LogLikelihoodRatio currResult = computeLikelihoodGivenCounts((double) (currAffectedCase+1), (double) (currAffectedControl+1), (double) (study.affectedSize+2), (double) (study.unaffectedSize+2), true);
				if(this.aboveThresholdA == -1 || (currResult.LLAlternative-currResult.LLNull) > (this.LLAlternative-this.LLNull)){
					this.aboveThresholdA = currAffectedCase;
					this.aboveThresholdU = currAffectedControl;
					this.cutoff = currCutOff;
					this.LLAlternative = currResult.LLAlternative;
					this.LLNull = currResult.LLNull;
				}
				currCutOff = score;
				if(study.sampleId2affected.get(sampleId)) currAffectedCase++;
				else currAffectedControl++;
			} else {
				System.err.println("BurdenLikelihoodRatio::computeLikelihoodRatioByBestCutOff -- SCHISM scores are not sorted in descending order.");
			}
		}
		
		//compute for the minimal burden cutoff
		LogLikelihoodRatio currResult = computeLikelihoodGivenCounts((double) (currAffectedCase+1), (double) (currAffectedControl+1), (double) (study.affectedSize+2), (double) (study.unaffectedSize+2), true);
		if(this.aboveThresholdA == -1 || (currResult.LLAlternative-currResult.LLNull) > (this.LLAlternative-this.LLNull)){
			this.aboveThresholdA = currAffectedCase;
			this.aboveThresholdU = currAffectedControl;
			this.cutoff = currCutOff;
			this.LLAlternative = currResult.LLAlternative;
			this.LLNull = currResult.LLNull;
		}
		
		//================ debug =================
		//System.out.println("best cutoff: "+String.valueOf(this.cutoff)+", CaseA: "+String.valueOf(this.aboveThresholdA)+", ControlA: "+String.valueOf(this.aboveThresholdU));
		//================ debug =================
		
	}
	
	//compute the Bernoulli log likelihood ratio between case-control only and both together
	LogLikelihoodRatio computeLikelihoodGivenCounts(double affectedCase, double affectedControl, double caseSize, double controlSize, boolean polarity){
		LogLikelihoodRatio result = new LogLikelihoodRatio();
		
		double allProb = ((double) (affectedCase+affectedControl))/((double) (caseSize+controlSize));
		
		double trueCaseProb = ((double) affectedCase)/((double) caseSize);
		double trueControlProb = ((double) affectedControl)/((double) controlSize);

		double caseProb = trueCaseProb;
		double controlProb = trueControlProb;

		//if the computation has direction (the default is "true" for burden statistic),
		//seeing a success in case always has a higher probability than that in control
		if(polarity){
			if(trueCaseProb < trueControlProb){
				caseProb = trueControlProb;
				controlProb = trueCaseProb;
			}
		}
		
		//here may have domain error exception for log calculation (log(0))
		//however, it is already dealt with by adding pseudo count when the function is called
		double caseLogLikelihood = affectedCase*Math.log(caseProb)+(caseSize-affectedCase)*Math.log(1-caseProb);
		double controlLogLikelihood = affectedControl*Math.log(controlProb)+(controlSize-affectedControl)*Math.log(1-controlProb);
		double allLogLikelihood = (affectedCase+affectedControl)*Math.log(allProb)+(caseSize+controlSize-affectedCase-affectedControl)*Math.log(1-allProb);

		result.LLAlternative = caseLogLikelihood+controlLogLikelihood;
		result.LLNull = allLogLikelihood;
		
		return result;
	}

}
