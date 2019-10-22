package Model;

import java.util.*;

import Data.*;

public class UnitStatistic {
	public PositionDistributionStatistic positionPart;
	public double burdenStatistic;
	public double positionStatistic;
	public double statistic;
	public int numOfExtreme;
	public int numOfExtremeBurden;
	public int numOfExtremePosition;
	
	public UnitStatistic(){
		this.burdenStatistic = 0.0;
		this.positionStatistic = 0.0;
		this.statistic = 0.0;
		this.numOfExtreme = 0;
		this.positionPart = null;
	}
	
	public UnitStatistic(CaseControlStudy study){
		this.burdenStatistic = 0.0;
		this.positionStatistic = 0.0;
		this.statistic = 0.0;
		this.numOfExtreme = 0;
		this.positionPart = null;
	}
	
	//compute burden statistic of the genomic unit for a case-control study
	public void computeBurdenStatistic(int unitId, CaseControlStudy study, Setting setting){
		//extract the list of variants associated with the tested genomic unit
		LinkedList<Integer> associatedVariantIds = study.genomicUnitId2genomicUnit.get(unitId).variantList;
		
		//compute individual burden
		Map<Integer, Double> sampleId2score = getAllIndividualBurden(associatedVariantIds,study.variantId2variant,setting.geneticModel);
		
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
	
	//compute position statistic of the genomic unit for a case-control study
	public void computePositionStatistic(int unitId, CaseControlStudy study, Setting setting){
		Iterator<Integer> variantIdItr = study.genomicUnitId2genomicUnit.get(unitId).variantList.iterator();
		LinkedList<Variant> variants = new LinkedList<Variant>();
		while(variantIdItr.hasNext()) variants.add(study.variantId2variant.get(variantIdItr.next()));
		
		LinkedList<Double> affectedCounts = new LinkedList<Double>();
		LinkedList<Double> unaffectedCounts = new LinkedList<Double>();
		
		Iterator<Variant> variantItr = variants.iterator();
		while(variantItr.hasNext()){
			Variant entryVar = variantItr.next();
			double affectedCount = 0;
			double unaffectedCount = 0;
			
			Iterator<Integer> sampleIdItr = entryVar.sampleIds.iterator();
			while(sampleIdItr.hasNext()){
				Sample entrySample = study.sampleId2sample.get(sampleIdItr.next());
				int genotype = entrySample.var2genotype.get(entryVar.id);
				double score = 0.0;
				
				switch(setting.geneticModel){
				case GeneticModel.additive:
					score = genotype*entryVar.preCompPositionScore;
					break;
				case GeneticModel.dominant:
					score = entryVar.preCompPositionScore;
					break;
				default:
					score = genotype*entryVar.preCompPositionScore;
					break;
				}
				
				if(study.sampleId2affected.get(entrySample.id)) affectedCount += score;
				else unaffectedCount += score;
			}
			
			affectedCounts.add(affectedCount);
			unaffectedCounts.add(unaffectedCount);
		}
		
		GenomicUnit entry = study.genomicUnitId2genomicUnit.get(unitId);

		PositionDistributionStatistic initial = new PositionDistributionStatistic(affectedCounts,unaffectedCounts,variants,entry,setting);
		this.positionStatistic = initial.LLAlternative-initial.LLNull;
		this.positionPart = initial;
	}

}
