package Data;

import java.util.*;

public class Sample {
	public int id;
	public int group;
	public double phenoTrait;
	public String name;
	public Map<Integer,Integer> var2genotype;
	
	public Sample(int id, double phenoTrait, String name){
		this.id = id;
		this.phenoTrait = phenoTrait;
		this.name = name;
		this.var2genotype = new HashMap<Integer,Integer>();
		this.group = 0;
	}
	
	public Sample(int id, double phenoTrait, String name, int group){
		this.id = id;
		this.phenoTrait = phenoTrait;
		this.name = name;
		this.var2genotype = new HashMap<Integer,Integer>();
		this.group = group;
	}
	
//	public void addGenotypes(String[] genotypes){
//		var2genotype = new HashMap<Integer,Integer>();
//		for(int i=1;i<genotypes.length+1;i++){
//			int genotype = Integer.valueOf(genotypes[i-1]);
//			if(genotype < 0){
//				System.err.println("Negative value in genotype file!");
//				System.exit(0);
//			}
//			if(genotype != 0) var2genotype.put(i, Integer.valueOf(genotypes[i-1]));
//		}
//	}
	
	public boolean isAffected(){
		return (phenoTrait > 0 ? true : false);
	}
}
