package Data;

import java.util.*;

public class Variant {
	public int id;
	public int genomicUnitId;
	public int subGenomicUnitId;
	public int pos;
	public double funcScore;
	public double popAlleleFreq;
	public double preCompBurdenScore;
	public double preCompPositionScore;
	public LinkedList<Integer> sampleIds;
	public LinkedList<Integer> genotypes;
	
	public Variant(int id, int genomicUnitId, int subGenomicUnitId, int pos, double funcScore){
		this.id = id;
		this.genomicUnitId = genomicUnitId;
		this.subGenomicUnitId = subGenomicUnitId;
		this.pos = pos;
		this.funcScore = funcScore;
		this.sampleIds = new LinkedList<Integer>();
		this.genotypes = new LinkedList<Integer>();
	}
	
	public void setPopAlleleFreq(double popAlleleFreq){
		this.popAlleleFreq = popAlleleFreq;
	}
	
	public void precomputeScore(Setting setting){
		this.preCompBurdenScore = setting.getScore(setting.burdenScoring, this.funcScore, this.popAlleleFreq);
		this.preCompPositionScore = setting.getScore(setting.positionScoring, this.funcScore, this.popAlleleFreq);
	}

}
