package Data;

import java.io.File;
import java.util.*;

public class ComplexDiseaseStudy {
	public Map<Integer,Sample> sampleId2sample;
	public Map<Integer,Variant> variantId2variant;
	public Map<Integer,GenomicUnit> genomicUnitId2genomicUnit;
	public Map<String,LinkedList<Integer>> setName2unitList;
	public Set<Integer> genomicUnitSet;
	public Set<Integer> filteredVarSet;
	
	public ComplexDiseaseStudy(String variantFile, String genomicUnitFile, Setting setting){
		this.sampleId2sample = new HashMap<Integer,Sample>();
		this.variantId2variant = new HashMap<Integer,Variant>();
		this.genomicUnitId2genomicUnit = new HashMap<Integer,GenomicUnit>();
		this.setName2unitList = new HashMap<String,LinkedList<Integer>>();
		this.genomicUnitSet = new HashSet<Integer>();
		this.filteredVarSet = new HashSet<Integer>();
		
		Map<String,Integer> genomicUnitName2id = new HashMap<String,Integer>();

		//read genomic unit (format: genomic_unit,length,subunit_list)
		try{
			Set<String> genomicUnitSet = new HashSet<String>();
			int genomicUnitId = 1;
			
			File fd = new File(genomicUnitFile);
			Scanner freader = new Scanner(fd);
			
			while(freader.hasNextLine()){
				String line = freader.nextLine();
				String[] fields = line.split(",");
				
				if(genomicUnitSet.contains(fields[0])){
					System.err.println("Error: duplicate genomic unit ("+fields[0]+")");
					System.exit(0);
				}
				
				genomicUnitSet.add(fields[0]);
				
				GenomicUnit one;
				try{
					if(fields.length == 2){
						one = new GenomicUnit(genomicUnitId,fields[0],Integer.valueOf(fields[1]),null);
					} else {
						if(fields[2].equals("")) one = new GenomicUnit(genomicUnitId,fields[0],Integer.valueOf(fields[1]),null);
						else {
							String[] subUnits = fields[2].split(";");
							one = new GenomicUnit(genomicUnitId,fields[0],Integer.valueOf(fields[1]),subUnits);
						}
					}
					this.genomicUnitId2genomicUnit.put(genomicUnitId, one);
					genomicUnitName2id.put(one.name, genomicUnitId);
				} catch (Exception e){
					System.err.println(e.toString());
					System.err.println("Genomic unit length is not a valid integer.");
					System.exit(0);
				}
				genomicUnitId++;
			}
			freader.close();
		} catch (Exception e){
			System.err.println(e.toString());
			System.err.println("Errors occur while reading genomic unit file (-u).");
			System.exit(0);
		}
		
		//read variants from file (format: variantName,genomic_unit,subunit,position,score)
		try{
			Set<String> variantSet = new HashSet<String>();
			int variantId = 1;
			
			File fd = new File(variantFile);
			Scanner freader = new Scanner(fd);
			
			while(freader.hasNextLine()){
				String line = freader.nextLine();
				String[] fields = line.split(",");
				
				if(variantSet.contains(fields[0])){
					System.err.println("Error: duplicate variant ("+fields[0]+")");
					System.err.println("Error occurs while reading variant file (-v).");
					System.exit(0);
				}
				
				variantSet.add(fields[0]);
				
				GenomicUnit oneGUnit = this.genomicUnitId2genomicUnit.get(genomicUnitName2id.get(fields[1]));
				if(oneGUnit == null){
					System.err.println("Genomic unit ("+fields[1]+") where variant ("+fields[0]+") locates does not exist.");
					System.err.println("Error occurs while reading variant file (-v).");
					System.exit(0);
				}
				
				//create this variant object
				Variant one = null;
				if(fields[2].equals(oneGUnit.name)){
					if(oneGUnit.length <= Integer.valueOf(fields[3])){
						System.err.println("WARNING: variant position ("+String.valueOf(Integer.valueOf(fields[3]))+") is outside the genomic unit "+oneGUnit.name+" (legnth="+String.valueOf(oneGUnit.length)+") so skipped ...");
						this.filteredVarSet.add(variantId);
//						System.exit(0);
					} else {
						//add this variant ID into variant list of the genomic unit which it belongs to
						oneGUnit.variantList.add(variantId);
					}
					one = new Variant(variantId,oneGUnit.id,-1,Integer.valueOf(fields[3]),Double.valueOf(fields[4]));
				} else {
					int i=0;
					for(i=0;i<oneGUnit.subUnits.length;i++){
						if(fields[2].equals(oneGUnit.subUnits[i])){
							one = new Variant(variantId,oneGUnit.id,i,Integer.valueOf(fields[3]),Double.valueOf(fields[4]));
							//add this variant ID into variant list of the genomic unit which it belongs to
							oneGUnit.variantList.add(variantId);
							break;
						}
					}
					if(i==oneGUnit.subUnits.length){
						System.err.println("WARNING: genomic unit ("+fields[1]+") does not have subunit ("+fields[2]+") in definition file so skipped ...");
						this.filteredVarSet.add(variantId);
//						one = new Variant(variantId,oneGUnit.id,oneGUnit.length,Integer.valueOf(fields[3]),Double.valueOf(fields[4]));
//						System.exit(0);
					}
				}
				
				this.variantId2variant.put(variantId, one);
				variantId++;
			}
			freader.close();
		} catch (Exception e){
			System.err.println(e.toString());
			System.err.println("Error occurs while reading variant file (-v).");
			System.exit(0);
		}
		
		//read sets of genomic units if specified
		if(setting.geneSetFile != null){
			try{
				File fd = new File(setting.geneSetFile);
				Scanner freader = new Scanner(fd);
				
				while(freader.hasNextLine()){
					String line = freader.nextLine();
					String[] fields = line.split("\t");
					
					if(fields.length < 3){
						System.err.println("WARNING: empty gene set ("+fields[0]+") so skipped ...");
						continue;
					}
					
					LinkedList<Integer> unitList = new LinkedList<Integer>();
					for(int i=2;i<fields.length;i++){
						if(genomicUnitName2id.containsKey(fields[i])){
							int genomicUnitId = genomicUnitName2id.get(fields[i]);
							
							//**********skip units with no variant
							if(this.genomicUnitId2genomicUnit.get(genomicUnitId).variantList.size() == 0) continue;
							
							unitList.add(genomicUnitId);
							this.genomicUnitSet.add(genomicUnitId);
						} else {
							System.err.println("WARNING: genomic unit, "+fields[i]+", is not defined so skipped ...");
						}
					}
					
					//*******modified
					this.setName2unitList.put(fields[0], unitList);
				}
				freader.close();
			} catch (Exception e){
				System.err.println(e.toString());
				System.err.println("Error occurs while reading gene set file (--geneSet).");
				System.exit(0);
			}
		}
	}

}
