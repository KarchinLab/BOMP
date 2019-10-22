package Data;

public class Setting {
	public String phenoFile;
	public String genoFile;
	public String variantFile;
	public String unitFile;
	public String outputFile;
	
	public int burdenScoring;
	public int positionScoring;
	//public double significanceLevel;
	public String geneSetFile;
	public boolean analysisByUnit;
	public boolean analysisBySet;
	public int studyType;
	public int geneticModel;
	public int burdenMethod;
	public int windowLength;
	public int windowOffset;
	public int windowGrowth;
	public boolean bestSegmentationByUniqueVar;
	public int numOfPermutation;
	
	public Setting(String[] args){
		this.phenoFile = null;
		this.genoFile = null;
		this.variantFile = null;
		this.unitFile = null;
		this.outputFile = "bomp.out";
		
		this.burdenScoring = 4;
		this.positionScoring = 1;
		//this.significanceLevel = 0.0;
		this.geneSetFile = null;	//null string means performing regular analysis for each genomic unit
		this.analysisByUnit = true;
		this.analysisBySet = false;
		this.studyType = ComplexDiseaseStudyType.DichotomousTrait;
		this.geneticModel = GeneticModel.additive;
		this.burdenMethod = StatisticalMethod.hardThresholding;
		this.windowLength = 8;
		this.windowOffset = 1;
		this.windowGrowth = 4;
		this.bestSegmentationByUniqueVar = true;
		this.numOfPermutation = 10000;
		
		if(args.length == 0){
			usage();
			System.exit(0);
		}
		
		try {
			for(int i=0;i<args.length;i++){
				if(args[i].equals("-p")){
					if(i == args.length-1){
						System.out.println("Missing argument after -p");
						usage();
						System.exit(0);
					} else {
						this.phenoFile = args[i+1];
						i++;
					}
				} else if(args[i].equals("-g")){
					if(i == args.length-1){
						System.out.println("Missing argument after -g");
						usage();
						System.exit(0);
					} else {
						this.genoFile = args[i+1];
						i++;
					}
				} else if(args[i].equals("-v")){
					if(i == args.length-1){
						System.out.println("Missing argument after -v");
						usage();
						System.exit(0);
					} else {
						this.variantFile = args[i+1];
						i++;
					}
				} else if(args[i].equals("-u")){
					if(i == args.length-1){
						System.out.println("Missing argument after -u");
						usage();
						System.exit(0);
					} else {
						this.unitFile = args[i+1];
						i++;
					}
				} else if(args[i].equals("-o")){
					if(i == args.length-1){
						System.out.println("Missing argument after -o");
						usage();
						System.exit(0);
					} else {
						this.outputFile = args[i+1];
						i++;
					}
				} else if(args[i].equals("--bScore")){
					if(i == args.length-1){
						System.out.println("Missing argument after --bScore");
						usage();
						System.exit(0);
					} else {
						this.burdenScoring = Integer.valueOf(args[i+1]);
						i++;
					}
				} else if(args[i].equals("--pScore")){
					if(i == args.length-1){
						System.out.println("Missing argument after --pScore");
						usage();
						System.exit(0);
					} else {
						this.positionScoring = Integer.valueOf(args[i+1]);
						i++;
					}
				} else if(args[i].equals("--geneSet")){
					if(i == args.length-1){
						System.out.println("Missing argument after --geneSet");
						usage();
						System.exit(0);
					} else {
						this.geneSetFile = args[i+1];
						this.analysisBySet = true;
						i++;
					}
				} else if(args[i].equals("--analyzeSingleUnit")){
					if(i == args.length-1){
						System.out.println("Missing argument after --analyzeSingleUnit");
						usage();
						System.exit(0);
					} else {
						if(args[i+1].toUpperCase().equals("F")) this.analysisByUnit = false;
						else this.analysisByUnit = true;
						i++;
					}
				} else if(args[i].equals("--geneticModel")){
					if(i == args.length-1){
						System.out.println("Missing argument after --geneticModel");
						usage();
						System.exit(0);
					} else {
						if(args[i+1].toUpperCase().equals("ADDITIVE")) this.geneticModel = GeneticModel.additive;
						else if(args[i+1].toUpperCase().equals("DOMINANT")) this.geneticModel = GeneticModel.dominant;
						i++;
					}
				} else if(args[i].equals("--burdenMethod")){
					if(i == args.length-1){
						System.out.println("Missing argument after --burdenMethod");
						usage();
						System.exit(0);
					} else {
						this.burdenMethod = Integer.valueOf(args[i+1]);
						i++;
					}
				} else if(args[i].equals("--windowSize")){
					if(i == args.length-1){
						System.out.println("Missing argument after --windowSize");
						usage();
						System.exit(0);
					} else {
						this.windowLength = Integer.valueOf(args[i+1]);
						i++;
					}
				} else if(args[i].equals("--segmentShift")){
					if(i == args.length-1){
						System.out.println("Missing argument after --segmentShift");
						usage();
						System.exit(0);
					} else {
						this.windowOffset = Integer.valueOf(args[i+1]);
						i++;
					}
				} else if(args[i].equals("--windowExpand")){
					if(i == args.length-1){
						System.out.println("Missing argument after --windowExpand");
						usage();
						System.exit(0);
					} else {
						this.windowGrowth = Integer.valueOf(args[i+1]);
						i++;
					}
				} else if(args[i].equals("--selectBestSegByAllVar")){
					this.bestSegmentationByUniqueVar = false;
				} else if(args[i].equals("--permutation")){
					if(i == args.length-1){
						System.out.println("Missing argument after --permutation");
						usage();
						System.exit(0);
					} else {
						this.numOfPermutation = Integer.valueOf(args[i+1]);
						i++;
					}
				} else {
					System.err.println("Unknown option: "+args[i]);
					System.exit(0);
				}
			}
		} catch (Exception e){
			System.err.println(e.toString());
			usage();
			System.exit(0);
		}
		
		if(this.phenoFile == null){
			System.err.println("Cannot find the phenotype file (-p).\n");
			usage();
			System.exit(0);
		}
		
		if(this.genoFile == null){
			System.err.println("Cannot find the genotype file (-g).\n");
			usage();
			System.exit(0);
		}
		
		if(this.variantFile == null){
			System.err.println("Cannot find the variant file (-v).\n");
			usage();
			System.exit(0);
		}
		
		if(this.unitFile == null){
			System.err.println("Cannot find the genomic unit file (-u).\n");
			usage();
			System.exit(0);
		}
		
		if(this.geneSetFile == null) this.analysisByUnit = true;
	}
	
	public void usage(){
		System.err.print("Usage:\n" +
				"java -jar bomp.jar -p <phenotype file> -g <genotype file> -v <variant file> -u <genomic unit file> <options>\n\n"+
				"Options:\n"+
				"-o <output file>\n" +
				"\tDefault output file is bomp.out\n\n"+
				"--bScore <integer number>\n" +
				"\t1:vallina score\n" +
				"\t2:functional score from variant file\n" +
				"\t3:frequency weight (1/sqrt(p(1-p)))\n" +
				"\t4:product of 2 and 3\n" +
				"\tDefault:4\n\n"+
				"--pScore <integer number>\n"+
				"\t1:vallina score\n" +
				"\t2:functional score from variant file\n" +
				"\t3:frequency weight (1/sqrt(p(1-p)))\n" +
				"\t4:product of 2 and 3\n" +
				"\tDefault:1\n\n"+
				"--geneSet <gene set file>\n" +
				"\tThe file where gene sets being analyzed are defined\n\n"+
				"--analyzeSingleUnit <T/F>\n" +
				"\tEnable(T)/Disable(F) analysis on every single genomic unit when a gene set\n" +
				"\tfile is specified.\n" +
				"\tIf the gene set file is not given, this option is always true.\n"+
				"\tDefault:T\n\n"+
				"--geneticModel <integer number>\n" +
				"\t1:additive model\n" +
				"\t2:dominant model\n" +
				"\tDefault:1\n\n"+
				"--burdenMethod <integer number>\n" +
				"\t1:hard thresholding by maximizing likelihood ratio\n" +
				"\tDefault:1\n\n"+
				"--windowSize <integer number>\n" +
				"\tThe initial window size of window segmentation used in position distribution\n" +
				"\tstatistic\n" +
				"\tDefault:8\n\n"+
				"--segmentShift <integer number>\n" +
				"\tThe amount of shift of window segmentation in position distribution statistic\n" +
				"\tThe amount of shift should be less than the initial window size.\n" +
				"\tDefault:1\n\n" +
				"--windowExpand <integer number>\n" +
				"\tHow many different window sizes will be tested to find the best window\n" +
				"\tsegmentation in position distribution statistic\n" +
				"\tThe window size and segmentation shift will be doubled each time.\n" +
				"\tDefault:4\n\n" +
				"--selectBestSegByAllVar\n" +
				"\tSelect the best window segmentation using all variants rather than unique\n" +
				"\tvariants in position distribution statistic calculation\n\n" +
				"--permutation <integer number>\n" +
				"\tNumber of permutations\n" +
				"\tDefault:10000\n"
				);
	}
	
	public double getScore(int mode, double funcScore, double popAlleleFreq){
		switch(mode){
		case 1: return 1.0;
		case 2: return funcScore;
		case 3: return 1.0/Math.sqrt(popAlleleFreq*(1-popAlleleFreq));
		case 4: return funcScore/Math.sqrt(popAlleleFreq*(1-popAlleleFreq));
		default: return 1.0;
		}
	}

}
