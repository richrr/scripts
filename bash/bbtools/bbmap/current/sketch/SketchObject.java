package sketch;

import java.util.Arrays;

import tax.TaxTree;

public class SketchObject {
	
	static synchronized void setTaxtree(String taxTreeFile){
		if(taxTreeFile==null){
			taxtree=null;
			return;
		}
		if(treefile!=null){
			assert(!treefile.equals(taxTreeFile));
			if(treefile.equals(taxTreeFile)){return;}
			treefile=taxTreeFile;
		}
		taxtree=TaxTree.loadTaxTree(taxTreeFile, System.err, hashNames);
	}
	
	static int parseMode(String arg, String a, String b){
		int mode_=-1;
		if(a.equals("mode")){
			if(b.equalsIgnoreCase("single") || b.equalsIgnoreCase("onesketch")){
				mode_=ONE_SKETCH;
			}else if(b.equalsIgnoreCase("taxa") || b.equalsIgnoreCase("pertaxa")){
				mode_=PER_TAXA;
			}else if(b.equalsIgnoreCase("sequence") || b.equalsIgnoreCase("persequence")){
				mode_=PER_SEQUENCE;
			}else if(b.equalsIgnoreCase("img")){
				mode_=IMG;
			}else if(Character.isDigit(b.charAt(0))){
				mode_=Integer.parseInt(b);
			}
		}else if(a.equalsIgnoreCase("single") || a.equalsIgnoreCase("onesketch")){
			mode_=ONE_SKETCH;
		}else if(a.equalsIgnoreCase("taxa") || a.equalsIgnoreCase("pertaxa")){
			mode_=PER_TAXA;
		}else if(a.equalsIgnoreCase("sequence") || a.equalsIgnoreCase("persequence")){
			mode_=PER_SEQUENCE;
		}
		return mode_;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Static Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public static final int RAW=0, HEX=1, A48=2;
	public static final char[] codingArray={'R', 'H', 'A'};
	
	public static int CODING=A48;
	public static boolean delta=true;
	public static float maxGenomeFraction=0.02f;
	public static boolean amino=false;
	
	static final byte[] hexTable=new byte[128];
	static {
		Arrays.fill(hexTable, (byte)-1);
		for(int i='0'; i<='9'; i++){
			hexTable[i]=(byte)(i-'0');
		}
		for(int i='A'; i<='F'; i++){
			hexTable[i]=hexTable[i+'a'-'A']=(byte)(i-'A'+10);
		}
		hexTable['x']=hexTable['X']=hexTable['-']=hexTable['+']=0;
	}
	
	public static TaxTree taxtree=null;
	private static String treefile=null;
	static boolean hashNames=false;
	static long maxReads=-1;
	
	public static final int ONE_SKETCH=1, PER_SEQUENCE=2, PER_TAXA=3, IMG=4;
	
}
