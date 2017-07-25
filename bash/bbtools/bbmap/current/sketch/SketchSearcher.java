package sketch;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Locale;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

import dna.Parser;
import kmer.AbstractKmerTable;
import kmer.KmerTableSet;
import kmer.HashBuffer;
import shared.Shared;
import shared.Tools;
import structures.Heap;
import structures.IntList;
import tax.TaxNode;
import tax.TaxTree;

public class SketchSearcher extends SketchObject {
	
	public SketchSearcher(){
		
	}

	public boolean parse(String arg, String a, String b){
		if(Parser.isJavaFlag(arg)){
			//do nothing
		}else if(a.equals("size") || a.equals("sketchsize")){
			sketchSize=Integer.parseInt(b);
		}else if(a.equals("verbose")){
			verbose=Tools.parseBoolean(b);
		}else if(a.equals("ref")){
			addFiles(b, refFiles);
		}else if(Sketch.parseCoding(a, b)){
			//Do nothing
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
		}else if(a.equals("threads") || a.equals("sketchthreads") || a.equals("t")){
			threads=Integer.parseInt(b);
		}
		
		else if(a.equals("maxfraction")){
			Sketch.maxGenomeFraction=Float.parseFloat(b);
		}else if(a.equals("k")){
			k=Integer.parseInt(b);
		}else if(a.equals("rcomp")){
			rcomp=Tools.parseBoolean(b);
		}
		
		else if(a.equals("level") || a.equals("taxlevel") || a.equals("minlevel")){
			taxLevel=Integer.parseInt(b);
		}else if(a.equals("printtax") || a.equals("printtaxa")){
			printTax=Tools.parseBoolean(b);
		}else if(a.equals("index") || a.equals("makeindex")){
			if(b!=null && "auto".equalsIgnoreCase(b)){
				autoIndex=true;
				makeIndex=true;
			}else{
				autoIndex=false;
				makeIndex=Tools.parseBoolean(b);
			}
		}else if(a.equals("indexsize") || a.equals("indexlimit")){
			indexLimit=Integer.parseInt(b);
		}
		
		else if(b==null && arg.indexOf('=')<0){//if(new File(arg).exists())
			refFiles.add(arg);
		}else{
			return false;
		}
		return true;
	}
	
	public void postParse(){
		if(amino){k=Tools.min(k, 12);}
	}

	public void compare(ArrayList<Sketch> sketches, StringBuilder sb){
		ConcurrentHashMap<Integer, Comparison> map=new ConcurrentHashMap<Integer, Comparison>();

		@SuppressWarnings("unchecked")
		ArrayList<Comparison>[] alca=new ArrayList[sketches.size()];
		
		final int[] buffer=Sketch.makeBuffer();
		AtomicInteger fakeID=new AtomicInteger(1900000000);
		for(int i=0; i<sketches.size(); i++){
			fakeID.set(1900000000);
			Sketch a=sketches.get(i);
			if(tables!=null){
//				System.err.println("foo: processing "+a.name()+", sketches.length="+sketches.size());
//				new Exception().printStackTrace();
				processUsingIndex(a, buffer, fakeID, map);
			}else if(threads<2){
				for(Sketch b : refSketches){
					processPair(a, b, buffer, fakeID, map);
				}
			}else{
				ArrayList<CompareThread> alct=new ArrayList<CompareThread>(threads);
				for(int t=0; t<threads; t++){
					alct.add(new CompareThread(a, t, threads, fakeID, map));
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
				}
				alct=null;
			}
			
			ArrayList<Comparison> al=alca[i]=new ArrayList<Comparison>(map.size());
			for(Entry<Integer, Comparison> e : map.entrySet()){
				al.add(e.getValue());
			}
			Shared.sort(al);
			Collections.reverse(al);
			while(al.size()>maxRecords){
				al.remove(al.size()-1);
			}
			map.clear();
		}
		
		for(int i=0; i<alca.length; i++){
			Sketch s=sketches.get(i);
			ArrayList<Comparison> al=alca[i];
			writeResults(al, s, sb);
		}
	}
	
	private class CompareThread extends Thread {
		
		CompareThread(Sketch a_, int pid_, int incr_, AtomicInteger fakeID_, ConcurrentHashMap<Integer, Comparison> map_){
			a=a_;
			pid=pid_;
			incr=incr_;
			fakeID=fakeID_;
			map=map_;
		}
		
		public void run(){
			if(tables!=null){
				processUsingIndex(a, buffer, fakeID, map);
			}else{
				for(int i=pid; i<refSketches.size(); i+=incr){
					Sketch b=refSketches.get(i);
					processPair(a, b, buffer, fakeID, map);
				}
			}
		}
		
		final AtomicInteger fakeID;
		final ConcurrentHashMap<Integer, Comparison> map;
		final int[] buffer=Sketch.makeBuffer();
		final int incr;
		final int pid;
		final Sketch a;
		
	}
	
	private ArrayList<Sketch> getSketchesViaIndex(Sketch a){
		final int[] singleton=new int[1];
		final IntList idList=new IntList(Tools.min(sketchSize, indexLimit, 1000));
		AbstractKmerTable[] tableArray=tables.tables();
		for(long key : a.array){
			AbstractKmerTable set=tableArray[(int)(key%WAYS)];
//			System.err.println(set.getValue(key));
			final int[] ids=set.getValues(key, singleton);
//			System.err.println(Arrays.toString(ids));
			if(ids!=null){
				for(int id : ids){
					if(id>=0){
						idList.add(id-1);
					}
				}
			}
		}
//		System.err.println("idList.size:"+idList.size);
		if(idList.size<1){return null;}
		idList.sort();
		ArrayList<Sketch> list=new ArrayList<Sketch>(Tools.min(8, idList.size));

		int last=-1;
		int count=0;
		for(int i=0; i<idList.size; i++){
			int id=idList.get(i);
			if(id==last){
//				System.err.println("A: "+last+", "+id+", "+count+", "+minCount);
				count++;
			}else{
//				System.err.println("B: "+last+", "+id+", "+count+", "+minCount);
				if(last>-1 && (count>=minCount)){
					list.add(refSketches.get(last));
				}
				last=id;
				count=0;
			}
		}
		if(last>-1 && (count>=minCount)){
			list.add(refSketches.get(last));
		}
//		idList.shrinkToUnique();
//		for(int i=0; i<idList.size; i++){
//			int id=idList.get(i);
//		}
		return list;
	}
	
	void processUsingIndex(Sketch a, int[] buffer, AtomicInteger fakeID, ConcurrentHashMap<Integer, Comparison> map){
//		Timer t=new Timer();
//		t.start("Began query.");
		ArrayList<Sketch> hits=getSketchesViaIndex(a);
//		t.stop("Got "+hits.size()+" hits.");
//		t.start();
//		System.err.println("hits: "+hits.size());
		if(hits==null || hits.isEmpty()){return;}
		comparisons.getAndAdd(hits.size());
		for(Sketch b : hits){
			processPair(a, b, buffer, fakeID, map);
		}
//		t.stop("Processed pairs.");
	}
	
	private void writeResults(ArrayList<Comparison> al, Sketch s, StringBuilder sb){
		sb.append("\nResults for "+s.name()+":\n\n");
		
		ArrayList<TaxNode> tnl=new ArrayList<TaxNode>();
		for(Comparison c : al){

			if(format==0){
				sb.append(String.format(Locale.ROOT, "WKID %.2f%%\tKID %.2f%%\tmatches %d\tcompared %d",
						100*c.idWeighted(), 100*c.idMin(), c.hits, c.minIndex)+"\ttaxID "+c.taxID+"\tgSize "+c.genomeSize+"\t"+c.name+"\n");
			}else{
				if(taxtree!=null && c.taxID>=0 && c.taxID<minFakeID){
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
				
				tnl.clear();
			}
		}
	}
	
	boolean processPair(Sketch a, Sketch b, int[] buffer, AtomicInteger fakeID, ConcurrentHashMap<Integer, Comparison> map){
//		System.err.println("Comparing "+a.name()+" and "+b.name());
		Comparison c=compareOneToOne(a, b, buffer, minCount, minID, null);
		if(c==null){return false;}
		if(c.taxID<0){c.taxID=fakeID.getAndIncrement();}
		
		TaxNode tn=(taxtree==null ? null : taxtree.getNode(b.taxID));
		if(tn!=null){
			c.name=tn.name;
			if(tn.level<taxLevel){
				tn=taxtree.getNodeAtLevel(b.taxID, taxLevel);
			}
		}
		Integer key=(tn==null ? c.taxID : tn.id);

		Comparison old=map.get(key);
//		System.err.println("A. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
		if(old!=null && old.compareTo(c)>0){return false;}
		
		old=map.put(key, c);
		while(old!=null && old.compareTo(c)>0){
//			System.err.println("B. Old: "+(old==null ? 0 : old.hits)+", new: "+c.hits);
			c=old;
			old=map.put(key, c);
		}
		return true;
	}
	
	private static Comparison compareOneToOne(final Sketch a, final Sketch b, int[] buffer, int minCount, float minID, Heap<Comparison> heap){
//		assert(heap!=null); //Optional, for testing.
		final int matches=a.countMatches(b, buffer);
		assert(matches==buffer[0]);
		final int div=buffer[1];
		
		if(matches<minCount || matches/(float)div<minID){
//			System.err.print(".");
			return null;
		}
		if(heap!=null && !heap.hasRoom() && heap.peek().hits>matches){return null;}
		
//		System.err.print("*");
		Comparison c=new Comparison(buffer, b);
		if(heap==null || heap.add(c)){return c;}
		return null;
	}
	
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
	
	public void loadReferences() {
		makeTool();
		refSketches=tool.loadSketches_MT(SketchObject.PER_TAXA, refFiles);
//		System.err.println("Sketches: "+refSketches.get(0).name());
		if(makeIndex){
			makeIndex();
		}
	}
	
	public void makeTool(){
		if(tool==null){tool=new SketchTool(sketchSize, k, minCount, rcomp);}
	}
	
	public ArrayList<Sketch> loadSketchesFromString(String sketchString){
		return tool.loadSketchesFromString(sketchString);
	}
	
	/*--------------------------------------------------------------*/
	
	
	public void makeIndex(){
		assert(tables==null);
		tables=new KmerTableSet(new String[] {"ways="+WAYS, "tabletype="+AbstractKmerTable.ARRAYH}, 20);
		tables.allocateTables();
		spawnIndexThreads();
	}
	
	
	/** Spawn index threads */
	private void spawnIndexThreads(){
		
		//Do anything necessary prior to processing
		
		//Determine how many threads may be used
		final int threads=Shared.threads();
		ArrayList<IndexThread> alht=new ArrayList<IndexThread>(threads);
		AtomicInteger ai=new AtomicInteger(0);
		AtomicLong al=new AtomicLong(0);
		for(int i=0; i<threads; i++){
			alht.add(new IndexThread(ai, al));
		}
		
		//Start the threads
		for(IndexThread pt : alht){
			pt.start();
		}
		
		//Wait for completion of all threads
		boolean success=true;
		long codesProcessed=0;
		for(IndexThread pt : alht){
			
			//Wait until this thread has terminated
			while(pt.getState()!=Thread.State.TERMINATED){
				try {
					//Attempt a join operation
					pt.join();
					synchronized(pt){
						codesProcessed+=pt.codesProcessedT;
					}
				} catch (InterruptedException e) {
					//Potentially handle this, if it is expected to occur
					e.printStackTrace();
				}
			}
			success&=pt.success;
		}
		
		//Track whether any threads failed
		if(!success){errorState=true;}
		
		System.err.println("Indexed "+al+" hashcodes."); //For some reason codesProcessed is wrong.
		
		//Do anything necessary after processing
//		System.gc();
	}
	
	/*--------------------------------------------------------------*/
	
	public class IndexThread extends Thread {
		
		public IndexThread(AtomicInteger nextIndex_, AtomicLong keyCount_){
			table=new HashBuffer(tables.tables(), 1000, 31, true, false);
			nextIndex=nextIndex_;
			keyCount=keyCount_;
		}
		
		@Override
		public void run(){
//			System.err.println("Thread running.");
			int id=nextIndex.getAndIncrement();
			final int numSketches=refSketches.size();
//			System.err.println("numSketches="+numSketches);
			while(id<numSketches){
				final Sketch sk=refSketches.get(id);
				final long[] array=sk.array;
				final int limit=Tools.min(array.length, sketchSize, indexLimit);
//				System.err.println("limit="+limit);
				for(int i=0; i<limit; i++){
					long key=array[i];
					table.set(key, id+1);
					codesProcessedT++;
				}
				id=nextIndex.getAndIncrement();
			}
			long temp=table.flush();
			
			synchronized(this){
				codesProcessedT+=0;
				success=true;
				keyCount.getAndAdd(codesProcessedT);
//				if(codesProcessedT>0){System.err.println(codesProcessedT);}
			}
		}
		
		AtomicInteger nextIndex;
		AtomicLong keyCount;
		long codesProcessedT=0;
		HashBuffer table;
		boolean success=false;
		
	}
	
	/*--------------------------------------------------------------*/
	
	public KmerTableSet tables;
	
	public int WAYS=31;
	public boolean autoIndex=true;
	public boolean makeIndex=true;
	public int sketchSize=10000;
	public int indexLimit=Integer.MAX_VALUE;
	public int k=31;
	public boolean rcomp=true;
	public int minFakeID=1900000000;
	public int maxRecords=100;
	public float minID=0.0002f;
	public int minCount=3;
	public int format=0;
	public int taxLevel=TaxTree.SPECIES;
	public boolean printTax=false;
	public SketchTool tool=null;
	public ArrayList<Sketch> refSketches=new ArrayList<Sketch>();
	public ArrayList<String> refFiles=new ArrayList<String>();
	public int threads=Shared.threads();
	boolean verbose;
	boolean errorState=false;
	AtomicLong comparisons=new AtomicLong(0);
	
}
