package sketch;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
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
import stream.ByteBuilder;
import structures.Heap;
import tax.PrintTaxonomy;
import tax.TaxNode;
import tax.TaxServer;
import tax.TaxTree;

/**
 * Compares one or more input sketches to a set of reference sketches.
 * 
 * @author Brian Bushnell
 * @date July 29, 2016
 *
 */
public class SendSketch extends SketchObject {
	

	
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
		SendSketch cs=new SendSketch(args);
		
		//Run the object
		cs.process(t);
	}
	
	/**
	 * Constructor.
	 * @param args Command line arguments
	 */
	public SendSketch(String[] args){
		
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
			}else if(a.equals("mincount") || a.equals("minhits")  || a.equals("hits") || a.equals("matches")){
				minCount=Integer.parseInt(b);
			}else if(a.equals("minid") || a.equals("id")){
				minID=Float.parseFloat(b);
				if(minID>1){minID/=100;}
			}else if(a.equals("records")){
				maxRecords=Integer.parseInt(b);
				assert(maxRecords>=1) : "Max records must be at least 1.";
			}else if(a.equals("format")){
				format=Integer.parseInt(b);
			}else if(a.equals("address")){
				address=b;
			}
			
			
			else if(a.equals("size")){
				size=(int)Tools.parseKMG(b);
			}else if(a.equals("maxfraction")){
				Sketch.maxGenomeFraction=Float.parseFloat(b);
			}else if(a.equals("k")){
				k=Integer.parseInt(b);
			}else if(a.equals("rcomp")){
				rcomp=Tools.parseBoolean(b);
			}else if(a.equals("reads")){
				maxReads=Long.parseLong(b);
			}
			
			else if(a.equals("level") || a.equals("taxlevel") || a.equals("minlevel")){
				taxLevel=Integer.parseInt(b);
			}else if(a.equals("printtax") || a.equals("printtaxa")){
				printTax=Tools.parseBoolean(b);
			}
			
			else if(parseMode(arg, a, b)>-1){
				mode_=parseMode(arg, a, b);
			}
			
			else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}
			
			else{
				outstream.println("Unknown parameter "+args[i]);
				assert(false) : "Unknown parameter "+args[i];
			}
		}
		mode=mode_;
		if(amino){k=Tools.min(k, 12);}
		
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
		
		//Adjust the number of threads for input file reading
		if(!ByteFile.FORCE_MODE_BF1 && !ByteFile.FORCE_MODE_BF2 && Shared.threads()>2){
			ByteFile.FORCE_MODE_BF2=true;
		}
		
		ffout=FileFormat.testOutput(out, FileFormat.TEXT, null, false, overwrite, append, false);
		
		//Ensure that no file was specified multiple times
		if(!Tools.testForDuplicateFiles(true, in.toArray(new String[0]))){
			throw new RuntimeException("\nSome file names were specified multiple times.\n");
		}
		
		tool=new SketchTool(size, k, minCount, rcomp);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Outer Methods        ----------------*/
	/*--------------------------------------------------------------*/
	
	public void process(Timer t){
		Timer ttotal=new Timer();
		
		t.start();
		inSketches=tool.loadSketches_MT(mode, in);
		final int numLoaded=(inSketches.size());
		t.stop();
		outstream.println("Loaded "+numLoaded+" sketches in \t"+t);
		t.start();
		
		if(ffout==null){return;}
		TextStreamWriter tsw=new TextStreamWriter(ffout);
		tsw.start();
		
		ByteBuilder bb=new ByteBuilder();
		
		int cntr=0;
		for(Sketch sk : inSketches){
			sk.toBytes(bb);
			cntr++;
			if(cntr>=100 || bb.length>500000){ //Don't allow too much data in a single transaction
				byte[] message=bb.toBytes();
				bb.clear();
				String result=sendAndReceive(message);
				tsw.print(result);
				cntr=0;
			}
		}
		
		if(bb.length>0){
			byte[] message=bb.toBytes();
			bb.clear();
			String result=sendAndReceive(message);
			tsw.println(result);
		}
		tsw.println();
		
//		System.err.println("sending "+bb.toString());
		
		tsw.poisonAndWait();
		errorState&=tsw.errorState;
		
		t.stop();
//		outstream.println("\nRan "+(inSketches.size()*refSketches.size())+" comparisons in \t"+t);
		ttotal.stop();
		outstream.println("Total Time: \t"+ttotal);
	}
	
	public String sendAndReceive(byte[] message){
		URL url=null;
		InputStream is=null;
		URLConnection connection=null;
		OutputStream os=null;
		try {
			url=new URL(address);
			connection=url.openConnection();
			connection.setDoOutput(true);
			os = connection.getOutputStream();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		try {
			os.write(message);
			is=connection.getInputStream();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String result=TaxServer.readStream(is);

		try {
			os.close();
			is.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return result;
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
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	private ArrayList<String> in=new ArrayList<String>();
	
	private String out="stdout.txt";

	private int size=10000;
	private int k=31;
	private boolean rcomp=true;
	private final SketchTool tool;
	private final int mode;
	
	private int maxRecords=100;
	
	private int minCount=3;
	
	private float minID=0.002f;
	
	private int taxLevel=2;
	
	private ArrayList<Sketch> inSketches;
	
	private boolean printTax=true;
	
	private int format=0;
	
	String address="https://taxonomy.jgi-psf.org/sketch";
	
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
