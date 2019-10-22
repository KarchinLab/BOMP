package Data;

import java.io.File;
import java.util.*;

import Model.FisherExact;

public class CaseControlStudy extends ComplexDiseaseStudy {
	public int affectedSize;
	public int unaffectedSize;
	public Map<Integer,Boolean> sampleId2affected;
	public Map<Integer,Integer> groupId2affectedSize;
	
	public CaseControlStudy(String phenoFile, String genoFile, String variantFile, String genomicUnitFile, Setting setting){
		super(variantFile,genomicUnitFile,setting);
		
		boolean hasGroup = false;
		
		this.affectedSize = 0;
		this.unaffectedSize = 0;
		this.sampleId2affected = new HashMap<Integer,Boolean>();
		this.groupId2affectedSize = null;
		
		//read phenotype file (format: sampleID,phenotype)
		try{
			Set<String> sampleSet = new HashSet<String>();
			int sampleId = 1;
			
			File fd = new File(phenoFile);
			Scanner freader = new Scanner(fd);
			
			while(freader.hasNextLine()){
				String line = freader.nextLine();
				String[] fields = line.split(",");
				
				if(sampleSet.contains(fields[0])){
					System.err.println("Error: duplicate sample ("+fields[0]+")");
					System.exit(0);
				}
				
				sampleSet.add(fields[0]);
				
				Sample one;
				if(fields.length < 3){
					if(hasGroup){
						System.err.println("Missing membership at line "+String.valueOf(sampleId));
						System.exit(0);
					}
					one = new Sample(sampleId,Double.valueOf(fields[1]),fields[0]);
				} else {
					if(sampleId > 1 && !hasGroup){
						System.err.println("Missing membership before line "+String.valueOf(sampleId));
						System.exit(0);
					}
					
					int group = Integer.valueOf(fields[2]);
					if(group == 0){
						System.err.println("sample group id cannot be 0.");
						System.exit(0);
					}
					one = new Sample(sampleId,Double.valueOf(fields[1]),fields[0],Integer.valueOf(fields[2]));
					hasGroup = true;
				}
				this.sampleId2sample.put(sampleId, one);
				sampleId++;
			}
			freader.close();
		} catch (Exception e){
			System.err.println(e.toString());
			System.err.println("Error occurs while reading phenotype file");
			System.exit(0);
		}
				
		//read genotype from file (format: genotype1,genotype2, ... for sample 1 in first row)
		try{
			File fd = new File(genoFile);
			Scanner freader = new Scanner(fd);
			
			int studySize = this.sampleId2sample.size();
			int variantSize = this.variantId2variant.size();
			int lineNo = 1;
			int sampleCount = 1;
			
			while(freader.hasNextLine()){
				String line = freader.nextLine();
				String[] fields = line.split(",");
				
				//check if the number of genotypes each row matches the number of variants defined in variant file
				if(fields.length != variantSize){
					System.err.println("Missing genotypes at line "+String.valueOf(lineNo));
					System.err.println("#variant in variant file: "+String.valueOf(variantSize)+", #variant in genotype file: "+String.valueOf(fields.length));
					System.err.println("Error occurs while reading genotype file (-g).");
					System.exit(0);
				}
				
				//build a variant-genotype map for each sample
				Sample entry = this.sampleId2sample.get(sampleCount);
				for(int i=0;i<fields.length;i++){
					int genotype = Integer.valueOf(fields[i]);
					if(genotype > 0 && !this.filteredVarSet.contains(i+1)){
						//variant-genotype map for sampleId = sampleCount
						entry.var2genotype.put(i+1, genotype);
						
						//build a list of samples who have the variant for each variant
						Variant varEntry = this.variantId2variant.get(i+1);
						varEntry.sampleIds.add(sampleCount);
						varEntry.genotypes.add(genotype);
					}
				}

				lineNo++;
				sampleCount++;
			}
			freader.close();
			
			//check if the number of samples (rows) matches the number of samples listed in phenotype file
			if(sampleCount-1 != studySize){
				System.err.println("Genotypes of some samples are missing!");
				System.err.println("#rows in phenotype file and genotype file don't match.");
				System.err.println("Error occurs while reading genotype file (-g).");
				System.exit(0);
			}
			
		} catch (Exception e){
			System.err.println(e.toString());
			System.err.println("Error occurs while reading genotype file (-g).");
			System.exit(0);
		}
		
		if(hasGroup) this.groupId2affectedSize = new HashMap<Integer,Integer>();
				
		//construct phenotype
		//for case-control study, 1.0 means affected while 0.0 means unaffected. otherwise, error occurs
		//count #affected samples in each group if membership is specified
		int sampleSize = this.sampleId2sample.size();
		for(int i=1;i<=sampleSize;i++){
			Sample one = this.sampleId2sample.get(i);
			
			if(hasGroup && !this.groupId2affectedSize.containsKey(one.group))
				this.groupId2affectedSize.put(one.group, 0);
			
			if(one.phenoTrait == 1.0){
				this.sampleId2affected.put(i, true);
				this.affectedSize++;
				if(hasGroup){
					int groupCount = this.groupId2affectedSize.get(one.group);
					this.groupId2affectedSize.put(one.group, groupCount+1);
				}
			} else if(one.phenoTrait == 0.0){
				this.sampleId2affected.put(i, false);
				this.unaffectedSize++;
			} else {
				System.err.println("Error: phenotypic status for case-control study must be 1 (case) or 0 (control)");
				System.exit(0);
			}
		}
		
		System.out.println("#genomic units: "+String.valueOf(this.genomicUnitId2genomicUnit.size()));
		System.out.println("#variants: "+String.valueOf(this.variantId2variant.size()));
		System.out.println("#affected: "+String.valueOf(this.affectedSize));
		System.out.println("#unaffected: "+String.valueOf(this.unaffectedSize));
		
		//estimate population allele frequency and precompute scores for each variant in case-control study
		this.preComputeVariatScore(setting);
	}
	
	public void preComputeVariatScore(Setting setting){
		//estimate population allele frequency and then precompute scores for each variant
		Iterator<Variant> varItr = this.variantId2variant.values().iterator();
		while(varItr.hasNext()){
			Variant varEntry = varItr.next();
			
			//count the number of affected samples and unaffected samples who have this variant
			Iterator<Integer> sampleIdItr = varEntry.sampleIds.iterator();
			Iterator<Integer> genotypeItr = varEntry.genotypes.iterator();
			int affectedAllele = 0;
			int unaffectedAllele = 0;
			while(sampleIdItr.hasNext()){
				int sampleId = sampleIdItr.next();
				int genotype = genotypeItr.next();
				if(this.sampleId2affected.get(sampleId)) affectedAllele += genotype;
				else unaffectedAllele += genotype;
			}
			
			//estimate population allele frequency for the variant (may have Type I error inflation so not used anymore)
//			if(unaffectedAllele >= affectedAllele){
//				//not a causal variant
//				varEntry.popAlleleFreq = ((double) (affectedAllele+unaffectedAllele+1))/((double) (2*(this.affectedSize+this.unaffectedSize)+2));
//			} else {
//				FisherExact pval = new FisherExact((this.affectedSize+this.unaffectedSize)*2);
//				double trust = pval.getTwoTailedP(affectedAllele, this.affectedSize*2-affectedAllele, unaffectedAllele, this.unaffectedSize*2-unaffectedAllele);
//				if(trust >= setting.significanceLevel){
//					//not a causal variant under desired significance level
//					varEntry.popAlleleFreq = ((double) (affectedAllele+unaffectedAllele+1))/((double) (2*(this.affectedSize+this.unaffectedSize)+2));
//				} else {
//					//it's a causal variant. estimate population allele frequency by posterior mean
//					double fromBoth = ((double) (affectedAllele+unaffectedAllele+1))/((double) (2*(this.affectedSize+this.unaffectedSize)+2));
//					double fromControlOnly = ((double) (unaffectedAllele+1))/((double) 2*this.unaffectedSize+2);
//					varEntry.popAlleleFreq = fromBoth*trust+fromControlOnly*(1-trust);
//				}
//			}
			
			double fromBoth = ((double) (affectedAllele+unaffectedAllele+1))/((double) (2*(this.affectedSize+this.unaffectedSize)+2));
//			double fromControlOnly = ((double) (unaffectedAllele+1))/((double) 2*this.unaffectedSize+2);
			
			varEntry.popAlleleFreq = fromBoth;
			
			//precompute scores
			varEntry.precomputeScore(setting);
		}
	}

}
