package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.Random;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

import dna.Parser;
import fileIO.ByteFile;
import fileIO.FileFormat;
import fileIO.ReadWrite;
import fileIO.ByteStreamWriter;
import kmer.HashArray1D;
import kmer.HashForest;
import kmer.KmerNode;
import kmer.KmerTableSet;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import stream.ByteBuilder;
import structures.LongHeap;
import structures.LongList;

/**
 * @author Brian Bushnell
 * @date June 28, 2016
 *
 */
public final class SketchTool extends SketchObject {
	
	/*--------------------------------------------------------------*/
	/*----------------         Main Method          ----------------*/
	/*--------------------------------------------------------------*/
	
	/**
	 * Code entrance from the command line.
	 * @param args Command line arguments
	 */
	public static void main(String[] args){
		
		args=Parser.parseConfig(args);
		if(Parser.parseHelp(args, true)){
			//printOptions();
			System.exit(0);
		}
		
		Timer t=new Timer();
		t.start();
		
		//Create a new CountKmersExact instance
		SketchTool mhs=new SketchTool(args);
		t.stop();
		System.err.println("Time: \t"+t);
	}
	
	public SketchTool(String[] args){
		System.err.println("Executing "+getClass().getName()+" "+Arrays.toString(args)+"\n");
		
		/* Set global defaults */
		ReadWrite.ZIPLEVEL=2;
		ReadWrite.USE_UNPIGZ=true;
		
		/* Initialize local variables with defaults */
		Parser parser=new Parser();
		
		ArrayList<String> list=new ArrayList<String>();
		
		int k_=31;
		int size_=10000;
		int mincount_=1;
		boolean rcomp_=true;
		float cutoff=0.02f;
		int mode=ONE_SKETCH;
		
		/* Parse arguments */
		for(int i=0; i<args.length; i++){

			final String arg=args[i];
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if("null".equalsIgnoreCase(b)){b=null;}
			while(a.charAt(0)=='-' && (a.indexOf('.')<0 || i>1 || !new File(a).exists())){a=a.substring(1);}
			
			if(a.equals("in")){
				if(b!=null){
					for(String s : b.split(",")){
						list.add(s);
					}
				}
			}else if(a.equals("k")){
				k_=Integer.parseInt(b);
			}else if(a.equals("size") || a.equals("length") || a.equals("len")){
				size_=Integer.parseInt(b);
			}else if(a.equals("mincount")){
				mincount_=Integer.parseInt(b);
			}else if(a.equals("cutoff")){
				cutoff=Float.parseFloat(b);
			}else if(a.equals("rcomp")){
				rcomp_=Tools.parseBoolean(b);
			}else if(a.equals("reads")){
				maxReads=Long.parseLong(b);
			}else if(Sketch.parseCoding(a, b)){
				//Do nothing
			}else if(parser.parse(arg, a, b)){
				//do nothing
			}else if(b==null){
				list.add(arg);
			}else{
				throw new RuntimeException("Unknown parameter "+args[i]);
			}
		}
		
		k=k_;
		size=size_;
		mincount=mincount_;
		rcomp=rcomp_;
		
		Timer t=new Timer();
		ArrayList<Sketch> sketches=loadSketches_MT(mode, list);
		t.stop();
		System.err.println("Loaded "+sketches.size()+" sketches in \t"+t);
		t.start();
		Sketch sketch=sketches.get(0);
		int[] buffer=Sketch.makeBuffer();
		for(int i=1; i<sketches.size(); i++){
			Sketch sketch2=sketches.get(i);
			final int matches=sketch.countMatches(sketch2, buffer);
			assert(matches==buffer[0]);
			final int len=buffer[1];
			final int min=buffer[2];
			final int max=buffer[3];
			
			if(matches>=mincount){
				final float idWeighted=matches/(float)len;
				final float idMin=matches/(float)min;
				final float idMax=matches/(float)max;
				
				if(idWeighted>=cutoff && matches>=mincount){
					System.out.println(String.format(Locale.ROOT, "%.2f%% WID, %.2f%% MINID, %.2f%% MAXID, %d matches, %d length",
							100*idWeighted, 100*idMin, 100*idMax, matches, len)+" for "+sketch.name()+" vs "+sketch2.name());
				}
			}
		}
		t.stop();
		System.err.println("Compared "+(sketches.size()-1)+" sketches in \t"+t);
		
//		String fname=list.get(0);
//		Sketch sketch=loadSketches(fname).get(0);
//		for(int i=1; i<list.size(); i++){
//			String fname2=list.get(i);
//			Sketch sketch2=loadSketches(fname2).get(0);
//			float identity=sketch.identity(sketch2);
//			System.out.println("Identity for "+fname+" vs "+fname2+":\t"+String.format(Locale.ROOT, "%.2f%%", 100*identity));
//		}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------      Normal Constructor      ----------------*/
	/*--------------------------------------------------------------*/
	
	public SketchTool(int size_, int k_, int mincount_, boolean rcomp_){
		k=k_;
		size=size_;
		mincount=mincount_;
		rcomp=rcomp_;

		assert(!amino || (k>0 && k<13)) : "Protein sketches require 0 < K < 13.";
		assert(k>0 && k<32) : "Sketches require 0 < K < 32.";
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/
	
	public Sketch toSketch(KmerTableSet tables, boolean multithreaded){
		final int threads=(multithreaded ? Tools.mid(1, Shared.threads(), tables.ways()) : 1);
		return (threads<2 ? toSketch_ST(tables) : toSketch_MT(tables, threads));
	}
	
	private Sketch toSketch_ST(KmerTableSet tables){
		SketchHeap heap=new SketchHeap(size);

		KmerTableSet kts=(KmerTableSet)tables;
		for(int tnum=0; tnum<kts.ways; tnum++){
			HashArray1D table=kts.getTable(tnum);
			toHeap(table, heap);
		}
		
		return new Sketch(heap);
	}
	
	private Sketch toSketch_MT(KmerTableSet tables, final int threads){
		ArrayList<SketchThread> alst=new ArrayList<SketchThread>(threads);
		AtomicInteger ai=new AtomicInteger(0);
		for(int i=0; i<threads; i++){
			alst.add(new SketchThread(ai, tables));
		}

		//Start the threads
		for(SketchThread pt : alst){
			pt.start();
		}

		ArrayList<SketchHeap> heaps=new ArrayList<SketchHeap>(threads);
		for(SketchThread pt : alst){

			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			if(pt.heap.size()>0){
				heaps.add(pt.heap);
			}
		}
		alst.clear();
		return toSketch(heaps);
	}
	
	public static final long[] toSketchArray(LongHeap heap, int maxLen){
		final int initial=heap.size();
		final int len=Tools.min(maxLen, initial);
		final long[] array=new long[len];
		
		int toSkip=heap.size()-len;
		for(int i=0; i<toSkip; i++){heap.poll();}
		for(int i=0; i<len; i++){
			array[i]=Long.MAX_VALUE-heap.poll();
		}
		Tools.reverseInPlace(array);
		assert(heap.size()==0 || initial>maxLen) : heap.size()+", "+len+", "+maxLen+", "+initial;
		heap.clear();
		return array;
	}
	
	public SketchHeap toHeap(HashArray1D table, SketchHeap heap){
//		if(heap==null){heap=new LongHeap(size, true);}
		long[] kmers=table.array();
		int[] counts=table.values();
		for(int i=0; i<table.arrayLength(); i++){
			int count=counts[i];
			if(count>=mincount){
				heap.genomeSize++;
				long hash=hash(kmers[i]);
				heap.add(hash);
			}
		}
		HashForest forest=table.victims();
		if(forest!=null){
			for(KmerNode kn : forest.array()){
				if(kn!=null){addRecursive(heap, kn);}
			}
		}
		return heap;
	}
	
//	public long[] toSketchArray(ArrayList<LongHeap> heaps){
//		if(heaps.size()==1){return toSketchArray(heaps.get(0));}
//		LongList list=new LongList(size);
//		for(LongHeap heap : heaps){
//			while(heap.size()>0){list.add(Long.MAX_VALUE-heap.poll());}
//		}
//		list.sort();
//		list.shrinkToUnique();
//		list.size=Tools.min(size, list.size);
//		return list.toArray();
//	}
	
	public Sketch toSketch(ArrayList<SketchHeap> heaps){
		SketchHeap a=heaps.get(0);
		for(int i=1; i<heaps.size(); i++){
			SketchHeap b=heaps.get(i);
			a.add(b);
		}
		return new Sketch(a);
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/
	
	private void addRecursive(SketchHeap heap, KmerNode kn){
		if(kn==null){return;}
		if(kn.count()>=mincount){
			heap.genomeSize++;
			long hash=hash(kn.pivot());
			heap.add(hash);
		}
		if(kn.left()!=null){addRecursive(heap, kn.left());}
		if(kn.right()!=null){addRecursive(heap, kn.right());}
	}
	
	/*--------------------------------------------------------------*/
	/*----------------             I/O              ----------------*/
	/*--------------------------------------------------------------*/
	
//	public static ArrayList<Sketch> loadSketches_ST(String...fnames){
//		ArrayList<Sketch> sketches=null;
//		for(String s : fnames){
//			ArrayList<Sketch> temp;
//			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
//				temp=loadSketches(s);
//			}else{
//				temp=loadSketches_ST(s.split(","));
//			}
//			if(sketches==null){sketches=temp;}
//			else{sketches.addAll(temp);}
//		}
//		return sketches;
//	}
	
//	public static ArrayList<Sketch> loadSketches_MT(ArrayList<String> fnames){
//		return loadSketches_MT(0, null, fnames.toArray(new String[0]));
//	}
	
	public ArrayList<Sketch> loadSketches_MT(int mode, ArrayList<String> fnames){
		return loadSketches_MT(mode, fnames.toArray(new String[0]));
	}
	
	public ArrayList<Sketch> loadSketches_MT(int mode, String...fnames){
		ConcurrentLinkedQueue<String> decomposedFnames=new ConcurrentLinkedQueue<String>();
		for(String s : fnames){
			if(s.indexOf(',')<0 || s.startsWith("stdin") || new File(s).exists()){
				decomposedFnames.add(s);
			}else{
				for(String s2 : s.split(",")){
					decomposedFnames.add(s2);
				}
			}
		}

		if(decomposedFnames.size()==0){return null;}
		if(decomposedFnames.size()==1){return loadSketches(decomposedFnames.poll(), null, mode);}
		
		
		//Determine how many threads may be used
		final int threads=Tools.min(Shared.threads(), decomposedFnames.size());

		//Fill a list with LoadThreads
		ArrayList<LoadThread> allt=new ArrayList<LoadThread>(threads);
		
		for(int i=0; i<threads; i++){
			allt.add(new LoadThread(decomposedFnames, mode));
		}
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		
		//Start the threads
		for(LoadThread lt : allt){lt.start();}

		//Wait for completion of all threads
		boolean success=true;
		for(LoadThread lt : allt){

			//Wait until this thread has terminated
			while(lt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					lt.join();
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			sketches.addAll(lt.list);
			success&=lt.success;
		}
		assert(success) : "Failure loading some files.";
		return sketches;
	}
	
	private class LoadThread extends Thread{
		
		public LoadThread(ConcurrentLinkedQueue<String> queue_, int mode_) {
			queue=queue_;
			list=new ArrayList<Sketch>();
//			map=map_;
//			if(map==null){list=new ArrayList<Sketch>();}
			smm=new SketchMakerMini(SketchTool.this, rcomp, mode_);
		}
		
		public void run(){
			success=false;
			for(String fname=queue.poll(); fname!=null; fname=queue.poll()){
				ArrayList<Sketch> temp=null;
				try {
					temp=loadSketches(fname, smm, smm.mode);
				} catch (Throwable e) {
					System.err.println("Failure loading "+fname+":\n"+e);
					success=false;
				}
				if(temp!=null){
					for(Sketch s : temp){add(s);}
				}
			}
			success=true;
		}
		
		private void add(Sketch s){
			if(list!=null){
				list.add(s);
				return;
			}
			assert(false) : "Unsupported."; //The map logic is broken; needs to be synchronized.
//			if(s.taxID<0){return;}
////			assert(s.taxID>-1) : s.toHeader();
//			TaxNode tn=tree.getNode(s.taxID);
//			while(tn!=null && tn.pid!=tn.id && tn.level<taxLevel){
//				TaxNode temp=tree.getNode(tn.pid);
//				if(temp==null){break;}
//				tn=temp;
//			}
//			if(tn==null){return;}
//			Integer key=tn.id;
//			Sketch old=map.get(key);
//			if(old==null){
//				s.taxID=key;
//				map.put(key, s);
//			}else{
//				synchronized(old){
//					old.add(s, maxLen);
//				}
//			}
		}
		
		final ConcurrentLinkedQueue<String> queue;
		ArrayList<Sketch> list;
		boolean success=false;
		final SketchMakerMini smm;
		
//		ConcurrentHashMap<Integer, Sketch> map;
		
	}
	
//	public ArrayList<Sketch> loadSketches(String fname){
//		
//	}
	
	public ArrayList<Sketch> loadSketches(String fname, SketchMakerMini smm, int mode){
		FileFormat ff=FileFormat.testInput(fname, FileFormat.TEXT, null, false, true);
		
		if(ff.fasta() || ff.fastq()){
			if(smm==null){smm=new SketchMakerMini(this, true, mode);}
			ArrayList<Sketch> sketches=smm.toSketches(fname);
			return sketches;
		}
		
		boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		ByteFile bf=ByteFile.makeByteFile(fname, false, false);
		int currentSketchSize=size;
		int taxID=-1;
		long imgID=-1;
		long genomeSize=0;
		String name=null, name0=null;
		LongList list=null;
		long sum=0;
		for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
			if(line.length>0){
//				System.err.println("Processing line "+new String(line));
				if(line[0]=='#'){
					if(list!=null){
						assert(list.size==list.array.length);
						list.shrink();
						Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSize, name, name0);
						sketches.add(sketch);
//						System.err.println("Made sketch "+sketch);
					}
					name=name0=null;
					list=null;
					sum=0;
					taxID=-1;
					imgID=-1;
					genomeSize=0;
					
					if(line.length>1){
						String[] split=new String(line, 1, line.length-1).split("\t");
						for(String s : split){
							if(s.startsWith("SIZE:")){//Sketch length
								currentSketchSize=Integer.parseInt(s.substring(5));
							}else if(s.startsWith("SZ:")){
								currentSketchSize=Integer.parseInt(s.substring(3));
							}else if(s.startsWith("CD:")){//Coding
								char c=s.charAt(3);
								A48=HEX=false;
								if(c=='A'){A48=true;}
								else if(c=='H'){HEX=true;}
								else{assert(c=='R');}
								delta=s.length()>4;
								assert(s.length()==4 || (s.length()==5 && s.charAt(4)=='D'));
							}else if(s.startsWith("GSIZE:")){//Genomic kmers
								genomeSize=Long.parseLong(s.substring(6));
							}else if(s.startsWith("GS:")){
								genomeSize=Long.parseLong(s.substring(3));
							}else if(s.startsWith("TAXID:")){
								taxID=Integer.parseInt(s.substring(6));
							}else if(s.startsWith("ID:")){
								taxID=Integer.parseInt(s.substring(3));
							}else if(s.startsWith("IMG:")){
								imgID=Long.parseLong(s.substring(4));
							}else if(s.startsWith("NAME:")){
								name=s.substring(5);
							}else if(s.startsWith("NM:")){
								name=s.substring(3);
							}else if(s.startsWith("NM0:")){
								name0=s.substring(4);
							}else{
								assert(false) : "Unsupported header tag "+s;
							}
						}
					}
					if(currentSketchSize>0){list=new LongList(currentSketchSize);}
				}else{
					long x=(A48 ? Sketch.parseA48(line) : HEX ? Sketch.parseHex(line) : Tools.parseLong(line));
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0) : x+"\n"+new String(line);
					assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
					list.add(delta ? sum : x);
					//						System.err.println("List="+list);
				}
			}
		}
		if(list==null){list=new LongList(0);}
		assert(list.size==list.array.length);
		list.shrink();
		Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSize, name, name0);
		sketches.add(sketch);
		return sketches;
	}
	
	public ArrayList<Sketch> loadSketchesFromString(String sketchString){
		boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
		
		ArrayList<Sketch> sketches=new ArrayList<Sketch>();
		int currentSketchSize=size;
		int taxID=-1;
		long imgID=-1;
		long genomeSize=0;
		String name=null, name0=null;
		LongList list=null;
		long sum=0;
		String[] split0=sketchString.split("\n");
		for(String line : split0){
			if(line.length()>0){
//				System.err.println("Processing line "+new String(line));
				if(line.charAt(0)=='#'){
					if(list!=null){
						assert(list.size==list.array.length);
						list.shrink();
						Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSize, name, name0);
						sketches.add(sketch);
//						System.err.println("Made sketch "+sketch);
					}
					name=name0=null;
					list=null;
					sum=0;
					taxID=-1;
					imgID=-1;
					genomeSize=0;
					
					if(line.length()>1){
						String[] split=line.substring(1).split("\t");
						for(String s : split){
							if(s.startsWith("SIZE:")){//Sketch length
								currentSketchSize=Integer.parseInt(s.substring(5));
							}else if(s.startsWith("SZ:")){
								currentSketchSize=Integer.parseInt(s.substring(3));
							}else if(s.startsWith("CD:")){//Coding
								char c=s.charAt(3);
								A48=HEX=false;
								if(c=='A'){A48=true;}
								else if(c=='H'){HEX=true;}
								else{assert(c=='R');}
								delta=s.length()>4;
								assert(s.length()==4 || (s.length()==5 && s.charAt(4)=='D'));
							}else if(s.startsWith("GSIZE:")){//Genomic kmers
								genomeSize=Long.parseLong(s.substring(6));
							}else if(s.startsWith("GS:")){
								genomeSize=Long.parseLong(s.substring(3));
							}else if(s.startsWith("TAXID:")){
								taxID=Integer.parseInt(s.substring(6));
							}else if(s.startsWith("ID:")){
								taxID=Integer.parseInt(s.substring(3));
							}else if(s.startsWith("IMG:")){
								imgID=Long.parseLong(s.substring(4));
							}else if(s.startsWith("NAME:")){
								name=s.substring(5);
							}else if(s.startsWith("NM:")){
								name=s.substring(3);
							}else if(s.startsWith("NM0:")){
								name0=s.substring(4);
							}else{
								assert(false) : "Unsupported header tag "+s;
							}
						}
					}
					if(currentSketchSize>0){list=new LongList(currentSketchSize);}
				}else{
					long x=(A48 ? Sketch.parseA48(line) : HEX ? Sketch.parseHex(line) : Long.parseLong(line));
//					System.err.println("sum="+sum+", x="+x+" -> "+(sum+x));
					sum+=x;
					assert(x>=0) : x+"\n"+new String(line);
					assert(sum>=0 || !delta) : "The sketch was made with delta compression off.  Please regenerate it.";
					list.add(delta ? sum : x);
					//						System.err.println("List="+list);
				}
			}
		}
		if(list==null){list=new LongList(1);}
		assert(list.size==list.array.length || list.size()==0);
		list.shrink();
		Sketch sketch=new Sketch(list.array, taxID, imgID, genomeSize, name, name0);
		sketches.add(sketch);
		return sketches;
	}
	
	public void write(ArrayList<Sketch> sketches, FileFormat ff[]){
		final int len=ff.length;
		ByteStreamWriter tsw[]=new ByteStreamWriter[len];
		for(int i=0; i<len; i++){
			tsw[i]=new ByteStreamWriter(ff[i]);
			tsw[i].start();
		}
		for(int i=0; i<sketches.size(); i++){
			write(sketches.get(i), tsw[i%len]);
		}
		for(int i=0; i<len; i++){
			tsw[i].poisonAndWait();
		}
	}
	
	public void write(ArrayList<Sketch> sketches, FileFormat ff){
		ByteStreamWriter tsw=new ByteStreamWriter(ff);
		tsw.start();
		for(Sketch sketch : sketches){
			write(sketch, tsw);
		}
		tsw.poisonAndWait();
	}
	
	public void write(Sketch sketch, FileFormat ff){
		ByteStreamWriter tsw=new ByteStreamWriter(ff);
		tsw.start();
		write(sketch, tsw);
		tsw.poisonAndWait();
	}
	
	public void write(Sketch sketch, ByteStreamWriter tsw){
		write(sketch, tsw, new byte[12]);
	}
	
	public void write(Sketch sketch, ByteStreamWriter tsw, byte[] temp){
		long prev=0;
		long[] array=sketch.array;
		final boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
		final ByteBuilder bb=sketch.toHeader().append('\n');
		for(int i=0; i<array.length; i++){
			long key=array[i];
			long x=key-prev;
			if(A48){
				Sketch.appendA48(x, bb, temp);
				bb.append('\n');
			}else if(HEX){
				bb.append(Long.toHexString(x));
				bb.append('\n');
			}else{
				bb.append(x);
				bb.append('\n');
			}
			if(delta){prev=key;}
		}
		tsw.print(bb);
	}
	
//	public void write(Sketch sketch, StringStreamWriter tsw, byte[] temp){
//		long prev=0;
//		long[] array=sketch.array;
//		final boolean A48=Sketch.CODING==Sketch.A48, HEX=Sketch.CODING==Sketch.HEX, delta=Sketch.delta;
////		final StringBuilder sb=sketch.toHeader().append('\n');
//		final StringBuilder bb=sketch.toHeader().append('\n');
//		tsw.print(sb.toString());
//		for(int i=0; i<array.length; i++){
//			long key=array[i];
//			long x=key-prev;
//			if(A48){
//				sb.setLength(0);
//				Sketch.appendA48(x, sb, temp);
//				sb.append('\n');
//				tsw.print(sb.toString());
//			}else if(HEX){
//				tsw.println(Long.toHexString(x));
//			}else{
//				sb.setLength(0);
//				sb.append(x);
//				sb.append('\n');
//				tsw.print(sb.toString());
//			}
//			if(delta){prev=key;}
//		}
//	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Nested Classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** Converts KmerTableSets to Heaps */
	private class SketchThread extends Thread {

		SketchThread(AtomicInteger next_, KmerTableSet kts_){
			next=next_;
			kts=kts_;
		}

		public void run(){
			final int ways=kts.ways();
			int tnum=next.getAndIncrement();
			while(tnum<ways){
				if(heap==null){heap=new SketchHeap(size);}
				HashArray1D table=kts.getTable(tnum);
				toHeap(table, heap);
				tnum=next.getAndIncrement();
			}
		}

		final AtomicInteger next;
		final KmerTableSet kts;
		SketchHeap heap;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Hashing            ----------------*/
	/*--------------------------------------------------------------*/
	
//	public long hashOld(long kmer){
//		assert(prime1>0);
//		long rot=Long.rotateRight(kmer^0x5555555555555555L, 9);
//		long a=rot%prime1;
//		long b=rot%prime2;
//		return Long.rotateRight(kmer, 17)^Long.rotateLeft(a, 31)^Long.rotateRight(b, 31);
//	}
	
	public static long[][] makeCodes(int symbols, int modes, long seed, boolean positive){
		Random randy=(seed>=0 ? new Random(seed) : new Random());
		long mask=positive ? Long.MAX_VALUE : -1L;
		long[][] r=new long[symbols][modes];
		for(int i=0; i<symbols; i++){
			for(int j=0; j<modes; j++){
				r[i][j]=randy.nextLong()&mask;
			}
		}
		return r;
	}
	
	public static final long hash(long kmer){
		long code=kmer;
//		assert(false) : Long.toBinaryString(middleMask);
		for(int i=0; i<8; i++){
			int x=(int)(kmer&0xFF);
			kmer>>=8;
			code^=codes[i][x];
//			code^=Long.rotateRight(codes[i][x], (int)(code&3)); //Might be better.
		}
		return code;
	}
	
	private static final long[][] codes=makeCodes(8, 256, 12345, false);
		
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	final int k;
	final int size;
	final boolean rcomp;
	final int mincount;
	
}
