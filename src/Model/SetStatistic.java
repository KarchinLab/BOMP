package Model;

import Data.*;

import java.util.*;

public class SetStatistic {
	public double burdenStatistic;
	public double positionStatistic;
	public double statistic;
	public int numOfExtreme;
	public int numOfExtremeBurden;
	public int numOfExtremePosition;
	
	public SetStatistic(){
		this.burdenStatistic = 0.0;
		this.positionStatistic = 0.0;
		this.statistic = 0.0;
		this.numOfExtreme = 0;
	}
	
	//compute burden set statistic and position set statistic at the same time by summing over
	//the statistic of each genomic unit in the set
	public void computeStatisticBySum(LinkedList<Integer> unitList, Map<Integer,UnitStatistic> unitStatMap){
		Iterator<Integer> unitIdItr = unitList.iterator();
		while(unitIdItr.hasNext()){
			UnitStatistic unitEntry = unitStatMap.get(unitIdItr.next());
			if(unitEntry != null){
				this.burdenStatistic += unitEntry.burdenStatistic;
				this.positionStatistic += unitEntry.positionStatistic;
				this.statistic += unitEntry.statistic;
			}
		}
	}
	
	//compute burden set statistic only by summing over the statistic of each genomic unit in the set
	public void computeBurdenStatistic(LinkedList<Integer> unitList, Map<Integer,UnitStatistic> unitStatMap){
		Iterator<Integer> unitIdItr = unitList.iterator();
		while(unitIdItr.hasNext()){
			UnitStatistic unitEntry = unitStatMap.get(unitIdItr.next());
			if(unitEntry != null){
				this.burdenStatistic += unitEntry.burdenStatistic;
			}
		}
	}
	
	//compute position set statistic only by summing over the statistic of each genomic unit in the set
	public void computePositionStatistic(LinkedList<Integer> unitList, Map<Integer,UnitStatistic> unitStatMap){
		double casesLogLikelihood = 0;
		double controlsLogLikelihood = 0;
		double samplesLogLikelihood = 0;
		double casesNormConstant = 0;
		double controlsNormConstant = 0;
		
		Iterator<Integer> unitIdItr = unitList.iterator();
		while(unitIdItr.hasNext()){
			UnitStatistic unitEntry = unitStatMap.get(unitIdItr.next());
			if(unitEntry != null){
				casesLogLikelihood += unitEntry.positionPart.caseUnormLogLikelihood;
				controlsLogLikelihood += unitEntry.positionPart.controlUnormLogLikelihood;
				samplesLogLikelihood += unitEntry.positionPart.allUnormLogLikelihood;
				casesNormConstant += unitEntry.positionPart.caseNormConstant;
				controlsNormConstant += unitEntry.positionPart.controlNormConstant;
				
				//this.positionStatistic += unitEntry.positionStatistic;
			}
		}
		
		if(casesNormConstant != 0) casesLogLikelihood -= casesNormConstant*Math.log(casesNormConstant);
		if(controlsNormConstant != 0) controlsLogLikelihood -= controlsNormConstant*Math.log(controlsNormConstant);
		if(casesNormConstant != 0 || controlsNormConstant != 0) samplesLogLikelihood -= (casesNormConstant+controlsNormConstant)*Math.log(casesNormConstant+controlsNormConstant);
		
		this.positionStatistic = casesLogLikelihood+controlsLogLikelihood-samplesLogLikelihood;
	}
	
	//compute burden set statistic for case control study
	public void computeBurdenStatistic(String setName, CaseControlStudy study, Setting setting){
		
		//put unique variants involved in the set of genomic units together
		Set<Integer> variantIdSet = new HashSet<Integer>();
		Iterator<Integer> genomicUnitItr = study.setName2unitList.get(setName).iterator();
		while(genomicUnitItr.hasNext()){
			Iterator<Integer> variantIdItr = study.genomicUnitId2genomicUnit.get(genomicUnitItr.next()).variantList.iterator();
			while(variantIdItr.hasNext()) variantIdSet.add(variantIdItr.next());
		}
		
		//compute individual burden
		Map<Integer,Double> sampleId2score = getAllIndividualBurden(new LinkedList<Integer>(variantIdSet),study.variantId2variant,setting.geneticModel);
		
		//compute burden statistic based on the specified statistical method
		switch(setting.burdenMethod){
		case StatisticalMethod.hardThresholding:
			BurdenHTStatistic one = new BurdenHTStatistic(sampleId2score,study);
			this.burdenStatistic = one.LLAlternative-one.LLNull;
			break;
		default:
			BurdenHTStatistic two = new BurdenHTStatistic(sampleId2score,study);
			this.burdenStatistic = two.LLAlternative-two.LLNull;
			break;
		}
	}
	
	//compute the burden score of the genomic unit for each individual
	Map<Integer,Double> getAllIndividualBurden(LinkedList<Integer> variantIds, Map<Integer,Variant> id2variant, int mode){
		//store individual burden in a map
		Map<Integer, Double> sampleId2burden = new TreeMap<Integer, Double>();
		
		//compute individual burden by aggregating all variants within the genomic unit
		Iterator<Integer> variantIdItr = variantIds.iterator();
		while(variantIdItr.hasNext()){
			int variantId = variantIdItr.next();
			Variant varEntry = id2variant.get(variantId);
			
			//aggregating the burden for each individual who has the variant
			Iterator<Integer> sampleIdItr = varEntry.sampleIds.iterator();
			Iterator<Integer> genotypeItr = varEntry.genotypes.iterator();
			while(sampleIdItr.hasNext()){
				int sampleId = sampleIdItr.next();
				int genotype = genotypeItr.next();
				double varBurden = 0.0;
				
				//compute burden based on genetic model
				switch(mode){
				case GeneticModel.additive:
					varBurden = genotype*varEntry.preCompBurdenScore;
					break;
				case GeneticModel.dominant:
					varBurden = varEntry.preCompBurdenScore;
					break;
				default:
					varBurden = genotype*varEntry.preCompBurdenScore;
					break;
				}
				
				//update individual burden
				if(sampleId2burden.containsKey(sampleId)){
					double score = sampleId2burden.get(sampleId);
					sampleId2burden.put(sampleId, score+varBurden);
				} else sampleId2burden.put(sampleId, varBurden);
			}
		}
		
		return sampleId2burden;
	}

}
