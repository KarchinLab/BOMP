package Data;

import Model.*;

public class CommandMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// TODO Auto-generated method stub

		String[] path = new String[5];
		path[0] = "/Users/yun-chingchen/Documents/projects/bomp/bompNew/bomp.new.phenotype";
		path[1] = "/Users/yun-chingchen/Documents/projects/bomp/bompNew/bomp.new.genotype";
		path[2] = "/Users/yun-chingchen/Documents/projects/bomp/bompNew/bomp.new.variant";
		path[3] = "/Users/yun-chingchen/Documents/projects/bomp/bompNew/bomp.new.unit";
		path[4] = "/Users/yun-chingchen/Documents/projects/bomp/bompNew/bomp.new.out";
		
		Setting setting = new Setting(args);
		Analyzer myAnalyzer = new Analyzer(setting,setting.phenoFile,setting.genoFile,setting.variantFile,setting.unitFile);
		myAnalyzer.output(setting.outputFile, setting);
//		Analyzer myAnalyzer = new Analyzer(setting,path[0],path[1],path[2],path[3]);
//		myAnalyzer.output(path[4], setting);
	}

}
