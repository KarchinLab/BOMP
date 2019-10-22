package Data;

import java.util.*;

public class GenomicUnit {
	public int id;
	public int length;
	public String name;
	public String[] subUnits;
	public LinkedList<Integer> variantList;
	
	public GenomicUnit(int id, String name, int length, String[] subUnits){
		this.id = id;
		this.length = length+1;	//because usually genomic or protein coordinate starting from 1
		this.name = name;
		this.subUnits = subUnits;
		this.variantList = new LinkedList<Integer>();
	}

}
