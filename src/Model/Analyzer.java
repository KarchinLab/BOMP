package Model;

import Data.*;

import java.io.FileWriter;
import java.util.*;

public class Analyzer {
	public Map<Integer,UnitStatistic> genomicUnitId2statistic;
	public Map<String,SetStatistic> set2statistic;
	public Map<Integer,UnitStatistic> genomicUnitId2permStatistic;
	public Set<Integer> targetedUnitSet;
	ComplexDiseaseStudy study;
	CaseControlStudy dichStudy;
	
	public Analyzer(Setting setting, String phenoFile, String genoFile, String variantFile, String genomicUnitFile){
		this.study = null;
		this.dichStudy = null;
		this.genomicUnitId2permStatistic = new HashMap<Integer,UnitStatistic>();
		this.genomicUnitId2statistic = new HashMap<Integer,UnitStatistic>();
		this.set2statistic = new HashMap<String,SetStatistic>();
		
		switch(setting.studyType){
		case ComplexDiseaseStudyType.DichotomousTrait:
			this.dichStudy = new CaseControlStudy(phenoFile, genoFile, variantFile, genomicUnitFile, setting);
			study = (ComplexDiseaseStudy) dichStudy;
			break;
		default:
			System.err.println("Unknown study type.");
			System.exit(0);
		}
				
		if(setting.analysisByUnit){
			//*****************remove units with no variant
			this.targetedUnitSet = new TreeSet<Integer>();
			
			Iterator<Integer> genomicUnitIdItr = study.genomicUnitId2genomicUnit.keySet().iterator();
			while(genomicUnitIdItr.hasNext()){
				int entry = genomicUnitIdItr.next();
				if(study.genomicUnitId2genomicUnit.get(entry).variantList.size() != 0)
					this.targetedUnitSet.add(entry);
				else {
					UnitStatistic one = new UnitStatistic();
					one.numOfExtreme = setting.numOfPermutation;
					one.numOfExtremeBurden = setting.numOfPermutation;
					one.numOfExtremePosition = setting.numOfPermutation;
					this.genomicUnitId2statistic.put(entry, one);
				}
			}
			
			if(setting.analysisBySet){
				LinkedList<String> removeList = new LinkedList<String>();
				Iterator<String> setNameItr = study.setName2unitList.keySet().iterator();
				while(setNameItr.hasNext()){
					String entry = setNameItr.next();
					if(study.setName2unitList.get(entry).size() == 0){
						SetStatistic one = new SetStatistic();
						one.numOfExtreme = setting.numOfPermutation;
						one.numOfExtremeBurden = setting.numOfPermutation;
						one.numOfExtremePosition = setting.numOfPermutation;
						this.set2statistic.put(entry, one);
						removeList.add(entry);
					}
				}
				setNameItr = removeList.iterator();
				while(setNameItr.hasNext()) study.setName2unitList.remove(setNameItr.next());
			}
			
			//************************
		} else if(setting.analysisBySet){
			this.targetedUnitSet = study.genomicUnitSet;
			
			//***********modified
			LinkedList<String> removeList = new LinkedList<String>();
			Iterator<String> setNameItr = study.setName2unitList.keySet().iterator();
			while(setNameItr.hasNext()){
				String entry = setNameItr.next();
				if(study.setName2unitList.get(entry).size() == 0){
					SetStatistic one = new SetStatistic();
					one.numOfExtreme = setting.numOfPermutation;
					one.numOfExtremeBurden = setting.numOfPermutation;
					one.numOfExtremePosition = setting.numOfPermutation;
					this.set2statistic.put(entry, one);
					removeList.add(entry);
				}
			}
			setNameItr = removeList.iterator();
			while(setNameItr.hasNext()) study.setName2unitList.remove(setNameItr.next());
			//***************************
		} else {
			System.err.println("Neither genomic unit analysis nor analysis of genomic unit sets. Nothing being analyzed.");
			System.exit(0);
		}
		
		//compute statistic for each genomic unit
		if(setting.analysisByUnit) this.computeUnitStatistic(setting);
		
		//compute set statistic if specified
		if(setting.analysisBySet) this.computeSetStatistic(setting);
		
		long startTime = System.nanoTime();
		long endTime = 0;
		
		//permutation
		List<Integer> sampleIdList = new ArrayList(study.sampleId2sample.keySet());
		for(int i=0;i<setting.numOfPermutation;i++){
			
			if(i%(setting.numOfPermutation/10) == 0){
				endTime = System.nanoTime();
				System.err.println(String.valueOf(i)+"th permutation");
				//System.err.println(String.valueOf(i)+"th permutation, duration: "+String.valueOf(endTime-startTime));
				//startTime = endTime;
			}
			
			//shuffle sample IDs
			//may shuffle within groups for population stratification
			Collections.shuffle(sampleIdList);
			Iterator<Integer> sampleIdItr = sampleIdList.iterator();
			
			//assign phenotype based on different types of studies
			switch(setting.studyType){
			case ComplexDiseaseStudyType.DichotomousTrait:
				if(dichStudy.groupId2affectedSize == null){
					//original permutation
					
					int index = 0;
					while(sampleIdItr.hasNext()){
						if(index < dichStudy.affectedSize) dichStudy.sampleId2affected.put(sampleIdItr.next(), true);
						else dichStudy.sampleId2affected.put(sampleIdItr.next(), false);
						index++;
					}
				} else {
					//stratified permutation if membership of stratified group is specified for each sample
					
					Map<Integer,Integer> groupId2counter = new HashMap<Integer,Integer>();
					Iterator<Integer> groupIdItr = dichStudy.groupId2affectedSize.keySet().iterator();
					while(groupIdItr.hasNext()) groupId2counter.put(groupIdItr.next(), 0);
					while(sampleIdItr.hasNext()){
						Sample one = dichStudy.sampleId2sample.get(sampleIdItr.next());
						int groupCounter = groupId2counter.get(one.group);
						int groupAffectedSize = dichStudy.groupId2affectedSize.get(one.group);
						if(groupCounter < groupAffectedSize){
							dichStudy.sampleId2affected.put(one.id, true);
							groupId2counter.put(one.group, groupCounter+1);
						} else dichStudy.sampleId2affected.put(one.id, false);
					}
				}
				
				//compute frequency weight again based on the current permutation
//				dichStudy.preComputeVariatScore(setting);
				break;
			default:
				break;
			}
			
			//compute statistic for each genomic unit
			if(setting.analysisByUnit) this.computePermutedUnitStatistic(setting);
			
			//compute set statistic if specified
			if(setting.analysisBySet) this.computePermutedSetStatistic(setting);
			
			//this.genomicUnitId2permStatistic = null;
			//this.genomicUnitId2permStatistic = new HashMap<Integer,UnitStatistic>();
		}
	}
	
	//compute statistic for each genomic unit individually
	void computeUnitStatistic(Setting setting){
		Iterator<Integer> genomicUnitItr = this.targetedUnitSet.iterator();
		while(genomicUnitItr.hasNext()){
			GenomicUnit unitEntry = study.genomicUnitId2genomicUnit.get(genomicUnitItr.next());
			UnitStatistic one = new UnitStatistic();
			switch(setting.studyType){
			case ComplexDiseaseStudyType.DichotomousTrait:
				one.computeBurdenStatistic(unitEntry.id, this.dichStudy, setting);
				one.computePositionStatistic(unitEntry.id, this.dichStudy, setting);
				one.statistic = one.burdenStatistic+one.positionStatistic;
				break;
			default:
				System.err.println("Unknown study type.");
				System.exit(0);
			}
			this.genomicUnitId2statistic.put(unitEntry.id, one);
		}
	}
	
	//compute set statistic
	//currently, the set burden statistic is calculated for the whole set
	//while the set position statistic is just the sum of individual position statistics in the set
	void computeSetStatistic(Setting setting){
		
		//if the individual position statistic is not calculated, then calculate
		if(!setting.analysisByUnit){
			Iterator<Integer> genomicUnitItr = this.targetedUnitSet.iterator();
			while(genomicUnitItr.hasNext()){
				GenomicUnit unitEntry = study.genomicUnitId2genomicUnit.get(genomicUnitItr.next());
				UnitStatistic one = new UnitStatistic();
				switch(setting.studyType){
				case ComplexDiseaseStudyType.DichotomousTrait:
					one.computePositionStatistic(unitEntry.id, this.dichStudy, setting);
					break;
				default:
					System.err.println("Unknown study type.");
					System.exit(0);
				}
				this.genomicUnitId2statistic.put(unitEntry.id, one);
			}
		}
		
		//compute set statistic for each set
		Iterator<String> setNameItr = this.study.setName2unitList.keySet().iterator();
		while(setNameItr.hasNext()){
			SetStatistic one = new SetStatistic();
			String setName = setNameItr.next();
			LinkedList<Integer> unitList = this.study.setName2unitList.get(setName);
			
			//compute set position statistic by summing over postion statistics of all individual genomic units in the set
			one.computePositionStatistic(unitList, this.genomicUnitId2statistic);
			
			//compute set burden statistic by considering the whole set of genomic units together
			switch(setting.studyType){
			case ComplexDiseaseStudyType.DichotomousTrait:
				one.computeBurdenStatistic(setName, this.dichStudy, setting);
				break;
			default:
				System.err.println("Unknown study type.");
				System.exit(0);
			}
			
			one.statistic = one.burdenStatistic+one.positionStatistic;
			this.set2statistic.put(setName, one);
		}
	}

	//compute statistic for each genomic unit individually
	void computePermutedUnitStatistic(Setting setting){
		Iterator<Integer> genomicUnitItr = this.targetedUnitSet.iterator();
		while(genomicUnitItr.hasNext()){
			GenomicUnit unitEntry = study.genomicUnitId2genomicUnit.get(genomicUnitItr.next());
			UnitStatistic one = new UnitStatistic();
			switch(setting.studyType){
			case ComplexDiseaseStudyType.DichotomousTrait:
				one.computeBurdenStatistic(unitEntry.id, this.dichStudy, setting);
				one.computePositionStatistic(unitEntry.id, this.dichStudy, setting);
				one.statistic = one.burdenStatistic+one.positionStatistic;
				break;
			default:
				System.err.println("Unknown study type.");
				System.exit(0);
			}
			
			this.genomicUnitId2permStatistic.put(unitEntry.id, one);
			UnitStatistic entryUnitStat = this.genomicUnitId2statistic.get(unitEntry.id);
			if(one.statistic >= entryUnitStat.statistic) entryUnitStat.numOfExtreme++;
			if(one.burdenStatistic >= entryUnitStat.burdenStatistic) entryUnitStat.numOfExtremeBurden++;
			if(one.positionStatistic >= entryUnitStat.positionStatistic) entryUnitStat.numOfExtremePosition++;
		}
	}
	
	//compute set statistic
	//currently, the set burden statistic is calculated for the whole set
	//while the set position statistic is just the sum of individual position statistics in the set
	void computePermutedSetStatistic(Setting setting){
		
		//if the individual position statistic is not calculated, then calculate
		if(!setting.analysisByUnit){
			Iterator<Integer> genomicUnitItr = this.targetedUnitSet.iterator();
			while(genomicUnitItr.hasNext()){
				GenomicUnit unitEntry = study.genomicUnitId2genomicUnit.get(genomicUnitItr.next());
				UnitStatistic one = new UnitStatistic();
				switch(setting.studyType){
				case ComplexDiseaseStudyType.DichotomousTrait:
					one.computePositionStatistic(unitEntry.id, this.dichStudy, setting);
					break;
				default:
					System.err.println("Unknown study type.");
					System.exit(0);
				}
				this.genomicUnitId2permStatistic.put(unitEntry.id, one);
			}
		}
		
		//compute set statistic for each set
		Iterator<String> setNameItr = this.study.setName2unitList.keySet().iterator();
		while(setNameItr.hasNext()){
			SetStatistic one = new SetStatistic();
			String setName = setNameItr.next();
			LinkedList<Integer> unitList = this.study.setName2unitList.get(setName);
			
			//compute set position statistic by summing over postion statistics of all individual genomic units in the set
			one.computePositionStatistic(unitList, this.genomicUnitId2permStatistic);
			
			//compute set burden statistic by considering the whole set of genomic units together
			switch(setting.studyType){
			case ComplexDiseaseStudyType.DichotomousTrait:
				one.computeBurdenStatistic(setName, this.dichStudy, setting);
				break;
			default:
				System.err.println("Unknown study type.");
				System.exit(0);
			}
			
			one.statistic = one.burdenStatistic+one.positionStatistic;
			SetStatistic entrySetStat = this.set2statistic.get(setName);
			if(one.statistic >= entrySetStat.statistic) entrySetStat.numOfExtreme++;
			if(one.burdenStatistic >= entrySetStat.burdenStatistic) entrySetStat.numOfExtremeBurden++;
			if(one.positionStatistic >= entrySetStat.positionStatistic) entrySetStat.numOfExtremePosition++;
		}
	}
	
	public void output(String outFileName, Setting setting){
		try{
			FileWriter outFile;
			
			if(setting.analysisByUnit){
				outFile = new FileWriter(outFileName+".gene");
				outFile.write("Name\t" +
						"statistic\t" +
						"p-value\t"+
						"burden statistic\t"+
						"bueden p-value\t"+
						"position statistic\t"+
						"position p-value\t"+
						"position best segment\n");
				
				Iterator<Integer> genomicUnitIdItr = this.genomicUnitId2statistic.keySet().iterator();
				while(genomicUnitIdItr.hasNext()){
					GenomicUnit entryUnit = study.genomicUnitId2genomicUnit.get(genomicUnitIdItr.next());
					UnitStatistic entryUnitStat = this.genomicUnitId2statistic.get(entryUnit.id);
					
					outFile.write(String.format("%s\t%.2f\t%.10f\t%.2f\t%.10f\t%.2f\t%.10f\t%d,%d\n", entryUnit.name,
							entryUnitStat.statistic,
							((double) (entryUnitStat.numOfExtreme+1))/((double) (setting.numOfPermutation+1)),
							entryUnitStat.burdenStatistic,
							((double) (entryUnitStat.numOfExtremeBurden+1))/((double) (setting.numOfPermutation+1)),
							entryUnitStat.positionStatistic,
							((double) (entryUnitStat.numOfExtremePosition+1))/((double) (setting.numOfPermutation+1)),
							entryUnitStat.positionPart.bestWindowSize,
							entryUnitStat.positionPart.bestOffset));
				}
				
				outFile.close();
			}
			
			
			if(setting.analysisBySet){
				outFile = new FileWriter(outFileName+".set");
				
				outFile.write("Name\t" +
						"statistic\t" +
						"p-value\t"+
						"burden statistic\t"+
						"bueden p-value\t"+
						"position statistic\t"+
						"position p-value\n");
				
				Iterator<String> setNameItr = this.set2statistic.keySet().iterator();
				while(setNameItr.hasNext()){
					String setName = setNameItr.next();
					SetStatistic entrySetStat = this.set2statistic.get(setName);
					
					outFile.write(String.format("%s\t%.2f\t%.10f\t%.2f\t%.10f\t%.2f\t%.10f\n", setName,
							entrySetStat.statistic,
							((double) (entrySetStat.numOfExtreme+1))/((double) (setting.numOfPermutation+1)),
							entrySetStat.burdenStatistic,
							((double) (entrySetStat.numOfExtremeBurden+1))/((double) (setting.numOfPermutation+1)),
							entrySetStat.positionStatistic,
							((double) (entrySetStat.numOfExtremePosition+1))/((double) (setting.numOfPermutation+1))));
				}
				
				outFile.close();
			}
			
			
		} catch (Exception e){
			System.err.println("Error occurs while writing results to file.");
			System.err.println(e);
		}
	}
}
