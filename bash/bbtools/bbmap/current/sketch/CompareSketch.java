package sketch;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

import dna.Parser;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import shared.ReadStats;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.Heap;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxTree;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class CompareSketch extends SketchObject {
	

	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		//Start a timer immediately upon code entrance.
		Timer t=new Timer();
		
		//Create an instance of this class
		CompareSketch cs=new CompareSketch(args);
		
		//Run the object
		cs.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public CompareSketch(String[] args){
		
		//Process any config files
		args=Parser.parseConfig(args);
		
		//Detect whether the uses needs help
		if(Parser.parseHelp(args, true)){
			printOptions();
			System.exit(0);
		}
		
		//Print the program name and arguments
		outstream.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		//Set some shared static variables regarding PIGZ
		ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=true;
		ReadWrite.MAX_ZIP_THREADS=Shared.threads();
		Shared.READ_BUFFER_LENGTH=1;
		
		int mode_=ONE_SKETCH;
		
		//Create a parser object
		Parser parser=new Parser();
		parser.out1="stdout.txt";
		
		//Parse each argument
		for(int i=0; i<args.length; i++){
			String arg=args[i];
			
			//Break arguments into their constituent parts, in the form of "a=b"
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			while(a.startsWith("-")){a=a.substring(1);} //Strip leading hyphens
			
			if(a.equals("verbose")){
				verbose=Tools.parseBoolean(b);
			}else if(a.equals("in")){
				addFiles(b, in);
			}else if(Sketch.parseCoding(a, b)){
				//Do nothing
			}else if(a.equals("parse_flag_goes_here")){
				long fake_variable=Tools.parseKMG(b);
				//Set a variable here
			}else if(a.equals("format")){
				format=Integer.parseInt(b);
			}else if(a.equals("reads")){
				maxReads=Long.parseLong(b);
			}
			
			else if(a.equals("taxtree") || a.equals("tree")){
				taxTreeFile=b;
				if("auto".equalsIgnoreCase(b)){taxTreeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("printtax") || a.equals("printtaxa")){
				printTax=Tools.parseBoolean(b);
			}
			
			else if(parseMode(arg, a, b)>-1){
				mode_=parseMode(arg, a, b);
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else if(searcher.parse(arg, a, b)){
				//do nothing
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		mode=mode_;
		searcher.postParse();
		
		{//Process parser fields
			overwrite=ReadStats.overwrite=parser.overwrite;
			append=ReadStats.append=parser.append;

			out=parser.out1;
		}
		
		//Ensure there is an input file
		if(in.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one input file is required.");
		}
		
		//Ensure there is an ref file
		if(searcher.refFiles.isEmpty()){
			printOptions();
			throw new RuntimeException("Error - at least one reference file is required.");
		}
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		
		//Ensure input files can be read
		if(!Tools.testInputFiles(false, true, taxTreeFile) || !Tools.testInputFiles(true, true, in.toArray(new String[0]))){
			throw new RuntimeException("\nCan't read to some input files.\n");
		}
		
		if(taxTreeFile!=null){setTaxtree(taxTreeFile);}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		
		outstream.println("Loading sketches.");
		searcher.makeTool();
		inSketches=searcher.tool.loadSketches_MT(mode, in);
		if(searcher.autoIndex){searcher.makeIndex=inSketches.size()>8;}
		searcher.loadReferences();
			
		final int numLoaded=(inSketches.size()+searcher.refSketches.size());
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketches in \t"+t);
		t.start();

		
		TextStreamWriter tsw=(ffout==null ? null : new TextStreamWriter(ffout));
		if(tsw!=null){tsw.start();}
		
		final int threads=Tools.min(Shared.threads(), inSketches.size());
		
		ArrayList<CompareThread> alct=new ArrayList<CompareThread>(threads);
		AtomicInteger next=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alct.add(new CompareThread(i, next, tsw));
		}
		for(CompareThread ct : alct){ct.start();}
		boolean success=true;
		for(CompareThread ct : alct){
			
			//Wait until this thread has terminated
			while(ct.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					ct.join();
				} catch (InterruptedException e) {
					e.printStackTrace();
				}
			}
			
			synchronized(ct){
				success&=ct.success;
			}
		}
		alct=null;
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		if(tsw!=null){tsw.poisonAndWait();}
		errorState&=tsw.errorState;
		
		t.stop();
		long comparisons=(searcher.makeIndex ? searcher.comparisons.get() : inSketches.size()*searcher.refSketches.size());
		outstream.println("\nRan "+comparisons+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	private void writeResults(ArrayList<Comparison> al, Sketch s, TextStreamWriter tsw){
		tsw.println("\nResults for "+s.name()+":\n");
		
		ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
		StringBuilder sb=new StringBuilder();
		for(Comparison c : al){

			if(format==0){
				tsw.print(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%\tmatches %d\tcompared %d",
						100*c.idWeighted(), 100*c.idMin(), c.hits, c.minIndex)+"\ttaxID "+c.taxID+"\tgSize "+c.genomeSize+"\t"+c.name+"\n");
				if(printTax && taxtree!=null){
					if(c.taxID>=0 && c.taxID<searcher.minFakeID){
						TaxNode tn=taxtree.getNode(c.taxID);
						if(tn!=null){
							PrintTaxonomy.printTaxonomy(tn, tsw, taxtree, TaxTree.DOMAIN);
						}
					}
					tsw.print("\n");
				}
			}else{
				if(taxtree!=null && c.taxID>=0 && c.taxID<searcher.minFakeID){
					TaxNode tn=taxtree.getNode(c.taxID);
					while(tn!=null && tn.pid!=tn.id && tn.level<=TaxTree.DOMAIN){
						tnl.add(tn);
						tn=taxtree.getNode(tn.pid);
					}
				}
				
				sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%\tmatches %d\tcompared %d\t",
						100*c.idWeighted(), 100*c.idMin(), c.hits, c.minIndex));
				sb.append("\ttaxID ").append(c.taxID).append('\t');
				sb.append(c.name).append('\t');
				
				for(int i=tnl.size()-1; i>=0; i--){
					TaxNode tn=tnl.get(i);
					sb.append(tn.name);
					if(i>0){sb.append(';');}
				}
				sb.append('\n');
				
				tsw.print(sb.toString());
				
				tnl.clear();
				sb.setLength(0);
			}
		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addFiles(String a, Collection<String> list){
		if(a==null){return;}
		File f=null;
		if(a.indexOf(',')>=0){f=new File(a);}
		if(f==null || f.exists()){
			list.add(a);
		}else{
			for(String s : a.split(",")){list.add(s);}
		}
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("Please run the shellscript with no arguments for usage information."); //TODO
	}
	
	public ArrayList<Comparison> toList(ConcurrentHashMap<Integer, Comparison> map){
		ArrayList<Comparison> al=new ArrayList<Comparison>(map.size());
		for(Entry<Integer, Comparison> e : map.entrySet()){
			al.add(e.getValue());
		}
		Shared.sort(al);
		Collections.reverse(al);
		while(al.size()>searcher.maxRecords){
			al.remove(al.size()-1);
		}
		return al;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Inner Classes        ----------------*/
	/*--------------------------------------------------------------*/
	
	private class CompareThread extends Thread {
		
		CompareThread(final int tid_, final AtomicInteger nextSketch_, TextStreamWriter tsw_){
			tid=tid_;
			nextSketch=nextSketch_;
			tsw=tsw_;
		}
		
		public void run(){
			success=false;
			final int inLim=inSketches.size();
			
			for(int inNum=nextSketch.getAndIncrement(); inNum<inLim; inNum=nextSketch.getAndIncrement()){
				Sketch a=inSketches.get(inNum);
				if(searcher.makeIndex){
					searcher.processUsingIndex(a, buffer, fakeID, map);
				}else{
					for(Sketch b : searcher.refSketches){
						searcher.processPair(a, b, buffer, fakeID, map);
					}
				}
				ArrayList<Comparison> al=toList(map);
				fakeID.set(searcher.minFakeID);
				map.clear();
				if(tsw!=null){
					synchronized(tsw){
						writeResults(al, a, tsw);
					}
				}
			}
			synchronized(this){success=true;}
		}
		
		private final int tid;
		private final int[] buffer=Sketch.makeBuffer();

		private final AtomicInteger nextSketch;
		private final AtomicInteger fakeID=new AtomicInteger(searcher.minFakeID);
		private ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>(101);
		final TextStreamWriter tsw;
		
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";
	
	private String taxTreeFile=null;

	private final int mode;
	
	private ArrayList<Sketch> inSketches;
	
	private boolean printTax=true;
	
	private int format=0;
	
	public final SketchSearcher searcher=new SketchSearcher();
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/

	/** Primary output file */
	private final FileFormat ffout;
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	/** Overwrite existing output files */
	private boolean overwrite=false;
	/** Append to existing output files */
	private boolean append=false;
	
}
