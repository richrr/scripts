package var2;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;

import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import shared.Shared;
import shared.Tools;
import stream.ByteBuilder;

public class VarMap {
	
	/*--------------------------------------------------------------*/
	/*----------------        Construction          ----------------*/
	/*--------------------------------------------------------------*/
	
	VarMap(ScafMap scafMap_){
		this(scafMap_, -1, -1, -1, -1, -1);
	}

	VarMap(ScafMap scafMap_, int ploidy_, float pairingRate_, float totalQualityAvg_, float mapqAvg_, float readLengthAvg_){
		scafMap=scafMap_;
		ploidy=ploidy_;
		properPairRate=pairingRate_;
		totalQualityAvg=totalQualityAvg_;
		totalMapqAvg=mapqAvg_;
		readLengthAvg=readLengthAvg_;
		maps=new ConcurrentHashMap[WAYS];
		for(int i=0; i<WAYS; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
	public static VarMap loadVars(String fname, ScafMap scafMap){
		final ByteFile bf=ByteFile.makeByteFile(fname, false, true);
		final VarMap varMap=new VarMap(scafMap);
		final byte delimiter='\t';
		int ploidy=-1;
		double pairingRate=-1;
		double mapqAvg=-1;
		double totalQualityAvg=-1;
		double readLengthAvg=-1;
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			if(line[0]!='#'){
				Var v=new Var(line, delimiter);
				varMap.addUnsynchronized(v);
			}else{
				String[] split=new String(line).split("\t");
				String a=split[0], b=(split.length>1 ? split[1] : null);
				assert(split.length>1) : new String(line);
				if(a.equalsIgnoreCase("#ploidy")){
					ploidy=Integer.parseInt(b);
				}else if(a.equalsIgnoreCase("#pairingRate")){
					pairingRate=Double.parseDouble(b);
				}else if(a.equalsIgnoreCase("#totalQualityAvg")){
					totalQualityAvg=Double.parseDouble(b);
				}else if(a.equalsIgnoreCase("#mapqAvg")){
					mapqAvg=Double.parseDouble(b);
				}else if(a.equalsIgnoreCase("#readLengthAvg")){
					readLengthAvg=Double.parseDouble(b);
				}
			}
			line=bf.nextLine();
		}
		bf.close();
		varMap.ploidy=ploidy;
		varMap.properPairRate=(float)pairingRate;
		varMap.totalQualityAvg=(float)totalQualityAvg;
		varMap.totalMapqAvg=(float)mapqAvg;
		varMap.readLengthAvg=(float)readLengthAvg;
		return varMap;
	}
	
	//#CHROM POS    ID        REF  ALT     QUAL
	public static VarMap loadVcf(String fname, ScafMap scafMap){
		ByteFile bf=ByteFile.makeByteFile(fname, false, true);
		VarMap varMap=new VarMap(scafMap);
		byte[] line=bf.nextLine();
		while(line!=null && line.length>0){
			if(line[0]!='#'){
				Var v;
				try {
					v = Var.fromVCF(line, scafMap);
					varMap.addUnsynchronized(v);
				} catch (Exception e) {
					System.err.println("Unable to parse VCF line: '"+new String(line)+"'");
//					throw new RuntimeException(e);
				}
			}else{
				String[] split=new String(line).split("=");
				if(split.length==2){
					String a=split[0], b=split[1];
					if(a.equalsIgnoreCase("##ploidy")){
						varMap.ploidy=Integer.parseInt(b);
					}else if(a.equalsIgnoreCase("##pairingRate")){
						varMap.properPairRate=(float) Double.parseDouble(b);
					}else if(a.equalsIgnoreCase("##totalQualityAvg")){
						varMap.totalQualityAvg=(float) Double.parseDouble(b);
					}else if(a.equalsIgnoreCase("##mapqAvg")){
						varMap.totalMapqAvg=(float) Double.parseDouble(b);
					}else if(a.equalsIgnoreCase("#readLengthAvg")){
						varMap.readLengthAvg=(float) Double.parseDouble(b);
					}
				}
			}
			line=bf.nextLine();
		}
		bf.close();
		return varMap;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	
	/*--------------------------------------------------------------*/
	/*----------------           Getters            ----------------*/
	/*--------------------------------------------------------------*/

	
	public boolean containsKey(Var v) {
		return get(v)!=null;
	}
	
	Var get(final Var v){
		final int way=v.start&MASK;
		return maps[way].get(v);
	}
	
	public long size(){
		long size=0;
		for(int i=0; i<maps.length; i++){size+=maps[i].size();}
		return size;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Adders            ----------------*/
	/*--------------------------------------------------------------*/
	
	private int add(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		synchronized(map){
			Var old=map.get(v);
			if(old==null){
				map.put(v, v);
				return 1;
			}
			else{
				synchronized(old){
					old.add(v);
				}
			}
		}
		return 0;
	}
	
	int addUnsynchronized(final Var v){
		final ConcurrentHashMap<Var, Var> map=maps[v.start&MASK];
		Var old=map.get(v);
		if(old==null){
			map.put(v, v);
			return 1;
		}
		old.add(v);
		return 0;
	}
	
	int dumpVars(HashMap<Var, Var> mapT){
		int added=0;
		@SuppressWarnings("unchecked")
		ArrayList<Var>[] absent=new ArrayList[WAYS];
		for(int i=0; i<WAYS; i++){
			absent[i]=new ArrayList<Var>();
		}
		for(Entry<Var, Var> e : mapT.entrySet()){
			Var v=e.getValue();
			final int way=v.start&MASK;
			ConcurrentHashMap<Var, Var> map=maps[way];
			Var old=map.get(v);
			if(old==null){absent[way].add(v);}
			else{
				synchronized(old){
					old.add(v);
				}
			}
		}
		
		mapT.clear();
		for(int way=0; way<WAYS; way++){
			ConcurrentHashMap<Var, Var> map=maps[way];
			ArrayList<Var> list=absent[way];
			synchronized(map){
				for(Var v : list){
					Var old=get(v);
					if(old==null){
						map.put(v, v);
						added++;
					}
					else{
						synchronized(old){
							old.add(v);
						}
					}
				}
			}
		}
		return added;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Other             ----------------*/
	/*--------------------------------------------------------------*/

	public long[] processVariantsST(VarFilter filter, long[] scoreArray) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		long[] types=new long[4];
		for(ConcurrentHashMap<Var, Var> map : maps){
			long[] types2=processVariants(map, filter, scoreArray, false);
			types2=processVariants(map, filter, scoreArray, true);
			types2=processVariants(map, filter, scoreArray, false);
			Tools.add(types, types2);
		}
		return types;
	}
	
	public long[] addSharedVariantsST(VarFilter filter, VarMap sharedVarMap) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		long[] types=new long[4];
		for(int i=0; i<maps.length; i++){
			long[] types2=addSharedVariants(maps[i], sharedVarMap.maps[i]);
			Tools.add(types, types2);
		}
		return types;
	}
	
	public long[] processVariantsMT(VarFilter filter, long[] scoreArray) {
		processVariantsMT_inner(filter, null, false);
		processVariantsMT_inner(filter, null, true);
		return processVariantsMT_inner(filter, scoreArray, false);
	}
	
	private long[] processVariantsMT_inner(VarFilter filter, long[] scoreArray, boolean processInsertions) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		ArrayList<ProcessThread> alpt=new ArrayList<ProcessThread>(WAYS);
		for(int i=0; i<WAYS; i++){
			ProcessThread pt=new ProcessThread(maps[i], filter, scoreArray!=null, processInsertions); 
			alpt.add(pt);
			pt.start();
		}
		
		long[] types=new long[4];
		boolean success=true;
		for(ProcessThread pt : alpt){
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
			
			//Accumulate per-thread statistics
			if(pt.types!=null){
				Tools.add(types, pt.types);
			}
			if(scoreArray!=null){Tools.add(scoreArray, pt.scoreArray);}
			success&=pt.success;
		}
		
		//Track whether any threads failed
//		if(!success){errorState=true;}
		
		return types;
	}
	
	private class ProcessThread extends Thread {
		
		ProcessThread(Map<Var, Var> map_, VarFilter filter_, boolean trackScores, boolean processInsertions_){
			map=map_;
			filter=filter_;
			scoreArray=(trackScores ? new long[200] : null);
			processInsertions=processInsertions_;
		}
		
		public void run(){
			types=processVariants(map, filter, scoreArray, processInsertions);
			success=true;
		}
		
		final VarFilter filter;
		final Map<Var, Var> map;
		long[] types;
		final long[] scoreArray;
		boolean processInsertions;
		boolean success=false;
	}

	private long[] processVariants(Map<Var, Var> map, VarFilter filter, long[] scoreArray, boolean processInsertions) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		Iterator<Entry<Var, Var>> iterator=map.entrySet().iterator();
		long[] types=new long[4];
		while(iterator.hasNext()){
			Entry<Var, Var> entry=iterator.next();
			final Var v=entry.getValue();
			
			if(processInsertions){
				assert(readLengthAvg>0);
				if(v.type()==Var.INS){
					synchronized(v){
						v.reviseAlleleFraction(readLengthAvg, scafMap.getScaffold(v.scafnum), this);
					}
				}
			}else{
				boolean pass=filter.passesFast(v);
				if(pass){
					v.calcCoverage(scafMap);
					pass=filter.passesFilter(v, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap);
				}
				if(pass){
					types[v.type()]++;
					if(scoreArray!=null){scoreArray[(int)v.phredScore(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap)]++;}
				}else{
					iterator.remove();
				}
			}
		}
		return types;
	}

	private long[] addSharedVariants(Map<Var, Var> map, Map<Var, Var> sharedMap) {
		assert(properPairRate>=0);
		assert(ploidy>0);
		assert(totalQualityAvg>=0);
		
		for(Var v : sharedMap.keySet()){
			if(!map.containsKey(v)){
				Var v2=new Var(v);
				map.put(v2, v2);
			}
		}
		
		long[] types=new long[4];
		for(Var v : sharedMap.keySet()){
			v.calcCoverage(scafMap);
			types[v.type()]++;
		}
		return types;
	}
	
	public Var[] toArray() {
		Var[] array=new Var[(int)size()];
		int i=0;
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Entry<Var, Var> e : map.entrySet()){
				array[i]=e.getValue();
				i++;
			}
		}
		return array;
	}
	
	public long[] countTypes() {
		long[] types=new long[4];
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Entry<Var, Var> e : map.entrySet()){
				types[e.getValue().type()]++;
			}
		}
		return types;
	}
	
	public void writeVarFile(FileFormat ff, VarFilter filter){
		Var[] array=toArray();
		Shared.sort(array);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder(33000);
		bb.append(Var.toHeader(properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.minAlleleFraction, filter.rarity, ploidy)).append('\n');
		for(Var v : array){
			v.toText(bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, filter.rarity, ploidy, scafMap);//TODO: Track depth
			bb.append('\n');
			if(bb.length()>16384){
				bsw.print(bb);
				bb.clear();
			}
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
		bsw.poisonAndWait();
	}
	
	public void writeVcfFile(String fname, VarFilter filter, String sampleName, long reads, long pairs, long properPairs, long bases, String ref, boolean trimWhitespace){
		FileFormat ff=FileFormat.testOutput(fname, FileFormat.TEXT, null, true, true, false, false);
		writeVcfFile(ff, filter, sampleName, reads, pairs, properPairs, bases, ref, trimWhitespace);
	}
	
	public void writeVcfFile(FileFormat ff, VarFilter filter, String sampleName, long reads, long pairs, long properPairs, long bases, String ref, boolean trimWhitespace){
		Var[] array=toArray();
		Shared.sort(array);
		ByteStreamWriter bsw=new ByteStreamWriter(ff);
		bsw.start();
		ByteBuilder bb=new ByteBuilder(33000);
		bb.append(Var.toVcfHeader(properPairRate, totalQualityAvg, totalMapqAvg, filter.rarity, filter.minAlleleFraction,
				ploidy, reads, pairs, properPairs, bases, scafMap, sampleName, ref, trimWhitespace)).append('\n');
		for(Var v : array){
			v.toVCF( bb, properPairRate, totalQualityAvg, totalMapqAvg, readLengthAvg, ploidy, scafMap, filter, trimWhitespace);
			bb.append('\n');
			if(bb.length()>16384){
				bsw.print(bb);
				bb.clear();
			}
		}
		if(bb.length()>0){
			bsw.print(bb);
			bb.clear();
		}
		bsw.poisonAndWait();
	}
	
	
	public void clear() {
		properPairRate=-1;
		pairedInSequencingRate=-1;
		totalQualityAvg=-1;
		totalMapqAvg=-1;
		readLengthAvg=-1;
		for(int i=0; i<maps.length; i++){
			maps[i]=new ConcurrentHashMap<Var, Var>();
		}
	}
	
	@Override
	public String toString(){
		StringBuilder sb=new StringBuilder();
		for(ConcurrentHashMap<Var, Var> map : maps){
			for(Var v : map.keySet()){
				sb.append(v.toString());
			}
		}
		return sb.toString();
	}
	
	/*--------------------------------------------------------------*/
	/*----------------           Fields             ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public int ploidy=-1;
	public float properPairRate=-1;
	public float pairedInSequencingRate=-1;
	public float totalQualityAvg=-1;
	public float totalMapqAvg=-1;
	public float readLengthAvg=-1;
	public final ScafMap scafMap;
	final ConcurrentHashMap<Var, Var>[] maps;
	
	/*--------------------------------------------------------------*/
	/*----------------        Static fields         ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	private static final int WAYS=4;
	public static final int MASK=WAYS-1;
	
}
