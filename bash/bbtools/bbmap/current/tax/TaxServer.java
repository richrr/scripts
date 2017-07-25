package tax;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.net.InetSocketAddress;

import com.sun.net.httpserver.HttpExchange;
import com.sun.net.httpserver.HttpHandler;
import com.sun.net.httpserver.HttpServer;

import dna.Parser;
import fileIO.ReadWrite;
import shared.Tools;
import sketch.Sketch;
import sketch.SketchObject;
import sketch.SketchSearcher;
import stream.KillSwitch;
import structures.IntList;

import com.sun.net.httpserver.Headers;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

/**
 * @author Shijie Yao, Brian Bushnell
 * @date Dec 13, 2016
 *
 */
public class TaxServer {

	public static void main(String[] args) throws Exception {
		TaxServer ts=new TaxServer(args);
		System.err.println("Ready!");
		//ts.begin();
	}
	
	public TaxServer(String[] args) throws Exception {
		int port_=3068;
		String killCode_=null;
		
		//Create a parser object
		Parser parser=new Parser();
		
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
			}else if(a.equals("verbose2")){
				verbose2=Tools.parseBoolean(b);
			}else if(a.equals("table") || a.equals("gi") || a.equals("gitable")){
				tableFile=b;
				if("auto".equalsIgnoreCase(b)){tableFile=TaxTree.defaultTableFile();}
			}else if(a.equals("tree") || a.equals("taxtree")){
				treeFile=b;
				if("auto".equalsIgnoreCase(b)){treeFile=TaxTree.defaultTreeFile();}
			}else if(a.equals("accession")){
				accessionFile=b;
				if("auto".equalsIgnoreCase(b)){accessionFile=TaxTree.defaultAccessionFile();}
			}else if(a.equals("reverse")){
				reverseOrder=Tools.parseBoolean(b);
			}else if(a.equals("domain")){
				domain=b;
				while(domain!=null && domain.endsWith("/")){domain=domain.substring(0, domain.length()-1);}
			}else if(a.equals("port")){
				port_=Integer.parseInt(b);
			}else if(a.equals("kill") || a.equals("killcode")){
				killCode_=b;
			}else if(searcher.parse(arg, a, b)){
				//do nothing
			}else if(parser.parse(arg, a, b)){//Parse standard flags in the parser
				//do nothing
			}else{
				throw new RuntimeException(arg);
			}
		}
		searcher.postParse();
		
		port=port_;
		killCode=killCode_;
		
		USAGE="Welcome to the JGI taxonomy server!\n"
				+ "This service provides taxonomy information from NCBI taxID numbers, gi numbers, organism names, and accessions.\n"
				+ "The output is formatted as a Json object.\n"
				+ "Usage:\n\n"
				+ "All addresses are assumed to be prefixed by "+domain+", e.g.\n"
				+ domain+"/tax/name/homo_sapiens\n"
				+ "\n"
				+ "/tax/name/homo_sapiens will give taxonomy information for an organism name.\n"
				+ "Names are case-insensitive and underscores are equivalent to spaces.\n"
				+ "/tax/taxid/9606 will give taxonomy information for an NCBI taxID.\n"
				+ "/tax/gi/1234 will give taxonomy information from an NCBI gi number.\n"
				+ "/tax/accession/NZ_AAAA01000057.1 will give taxonomy information from an accession.\n"
				+ "/tax/header/ will accept an NCBI sequence header such as gi|7|emb|X51700.1| Bos taurus\n"
				+ "Vertical bars (|) cause problems for curl and can be replaced by tilde (~).\n"
				+ "\nComma-delimited lists are accepted for bulk queries, such as tax/gi/1234,7000,42\n"
				+ "For plaintext (non-Json) results, use the pt_ or sc_ prefix.\n"
				+ "pt_ will give just the taxID, while sc_ will give the whole lineage, semicolon-delimited.\n"
				+ "For example:\n\n"
				+ "/tax/pt_name/homo_sapiens\n"
				+ "/tax/sc_gi/1234\n"
				+ "\nTo find the common ancestor of multiple organisms, add /ancestor/. For example:\n"
				+ "/tax/taxid/ancestor/1234,5678,42\n"
				+ "/tax/name/ancestor/homo_sapiens,canis_lupus,bos_taurus\n"
				+ "\nFor a simplified taxonomic tree, use /simpletax or /stax instead of /tax.\n"
				+ "This will ignore unranked or uncommon levels like tribe and parvorder, and only display the following levels:\n"
				+ "SUBSPECIES, SPECIES, GENUS, FAMILY, ORDER, CLASS, PHYLUM, KINGDOM, SUPERKINGDOM, DOMAIN\n"
				+ "For example:\n"
				+ "/simpletax/taxid/ancestor/1234\n"
				+ "\nTo print taxonomy from the command line in Linux, use curl:\n"
				+ "curl http://taxonomy.jgi-psf.org/tax/taxid/9606\n";
		
		typeMap=makeTypeMap();
		commonMap=makeCommonMap();
		
		
		if(tableFile!=null){
			outstream.println("Loading gi table.");
			GiToNcbi.initialize(tableFile);
		}
		if(treeFile!=null){
			outstream.println("Loading tree.");
			tree=ReadWrite.read(TaxTree.class, treeFile, true);
			if(tree.nameMap==null){
				outstream.println("Hashing names.");
				tree.hashNames();
			}
			assert(tree.nameMap!=null);
		}else{
			tree=null;
			throw new RuntimeException("No tree specified.");
		}
		if(accessionFile!=null){
			AccessionToTaxid.tree=tree;
			outstream.println("Loading accession table.");
			AccessionToTaxid.load(accessionFile);
			System.gc();
		}
		
//		searcher.threads=Shared.threads();
		SketchObject.taxtree=tree;
		
		if(!searcher.refFiles.isEmpty()){
			outstream.println("Loading sketches.");
			searcher.loadReferences();
			System.gc();
		}
		
		try {
			httpServer = HttpServer.create(new InetSocketAddress(port), 0);
		} catch (IOException e) {
			throw(e);
		}
		httpServer.createContext("/tax", new TaxHandler(false));
		httpServer.createContext("/stax", new TaxHandler(true));
		httpServer.createContext("/simpletax", new TaxHandler(true));
		httpServer.createContext("/sketch", new SketchHandler());
		if(killCode!=null){
			httpServer.createContext("/kill", new KillHandler());
		}
//		httpServer.createContext("/help", new HelpHandler());
		httpServer.createContext("/", new HelpHandler());
		httpServer.setExecutor(null); // creates a default executor
		httpServer.start();
	}

	class HelpHandler implements HttpHandler {
		
		public void handle(HttpExchange t) throws IOException {
			
			{
				Headers h = t.getResponseHeaders();
				String type="text/plain";
				h.add("Content-Type", type);
			}
			
			String response=USAGE;
			if(verbose2){System.out.println("Sending: "+response);}
			
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.getBytes());
			os.close();
		}
	}

	class KillHandler implements HttpHandler {
		
		public void handle(HttpExchange t) throws IOException {
			
			String rparam = t.getRequestURI().toString();
			while(rparam.startsWith("/")){
				rparam = rparam.substring(1);
			}
			while(rparam.endsWith("/")){
				rparam = rparam.substring(0, rparam.length()-1);
			}
			if(verbose){System.out.println(rparam);}

			String[] params = rparam.split("/");
			if(verbose2){System.out.println(Arrays.toString(params));}

			if(params.length>1){
				if(params[1].equals(killCode)){
					InetSocketAddress remote=t.getRemoteAddress();
					System.out.println("Killed by remote address "+remote);
					KillSwitch.killSilent();
				}
			}
			
			{
				Headers h = t.getResponseHeaders();
				String type="text/plain";
				h.add("Content-Type", type);
			}
			
			String response=BAD_CODE;
			if(verbose2){System.out.println("Sending: "+response);}
			
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.getBytes());
			os.close();
		}
	}

	class SketchHandler implements HttpHandler {
		
		public SketchHandler(){
		}
		
		public void handle(HttpExchange t) throws IOException {
			
			//String query = t.getRequestURI().getQuery(); //the KEY=VAL&KEY=VAL params in URL
			String rparam = t.getRequestURI().toString();   //restful style params, KEY/VAL in URL 
			while(rparam.startsWith("/")){
				rparam = rparam.substring(1);
			}
			while(rparam.endsWith("/")){
				rparam = rparam.substring(0, rparam.length()-1);
			}
			if(verbose){System.out.println(rparam);}
			
			assert(rparam.startsWith("sketch"));
			if(rparam.startsWith("sketch/")){rparam=rparam.substring(7);}
			else{rparam=rparam.substring(6);}

			if(verbose2){System.out.println(rparam);}
			
			boolean fileMode=false;
			if(rparam.startsWith("file/")){
				rparam=rparam.substring(5);
				fileMode=true;
			}

			if(verbose2){System.out.println(rparam);}
			if(verbose2){System.out.println("fileMode="+fileMode);}
			
			ArrayList<Sketch> sketches=null;
			if(fileMode){
				sketches=searcher.tool.loadSketches(rparam, null, SketchObject.ONE_SKETCH);
			}else{
				InputStream is=t.getRequestBody();
				String s=readStream(is);
//				System.err.println("Received "+s);
				if(s!=null && s.length()>0){
					sketches=searcher.loadSketchesFromString(s);
				}
			}
			
			StringBuilder response=new StringBuilder();
			if(sketches==null || sketches.isEmpty()){
				response.append("Error.");
			}else{
//				System.err.println("Received "+sketches.get(0).name()+", size "+sketches.get(0).array.length);
				searcher.compare(sketches, response);
//				System.err.println("Result: '"+response+"'");
			}
			
			{
				Headers h = t.getResponseHeaders();
				String type="text/plain";
				h.add("Content-Type", type);
			}
			
			if(verbose2){System.out.println("Sending: "+response);}
			
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.toString().getBytes());
			os.close();
		}
	}
	
	public static String readStream(InputStream is){
		try {
			byte[] buffer=new byte[256];
			int count=is.read(buffer);
			int next=0;
			while(count>-1){
				next+=count;
				if(next>=buffer.length){
					buffer=Arrays.copyOf(buffer, buffer.length*2);
				}
				count=is.read(buffer, next, buffer.length-next);
			}
			is.close();
			return new String(buffer, 0, next);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	class TaxHandler implements HttpHandler {
		
		public TaxHandler(boolean skipNonCanonical_){
			skipNonCanonical=skipNonCanonical_;
		}
		
		public void handle(HttpExchange t) throws IOException {
			
			//String query = t.getRequestURI().getQuery(); //the KEY=VAL&KEY=VAL params in URL
			String rparam = t.getRequestURI().toString();   //restful style params, KEY/VAL in URL 
			while(rparam.startsWith("/")){
				rparam = rparam.substring(1);
			}
			while(rparam.endsWith("/")){
				rparam = rparam.substring(0, rparam.length()-1);
			}
			if(verbose){System.out.println(rparam);}

			String[] params = rparam.split("/");
			if(verbose2){System.out.println(Arrays.toString(params));}
			
			final String response=toResponse(skipNonCanonical, params);
			
			{
				Headers h = t.getResponseHeaders();
				String type=response.startsWith("{") ? "application/json" : "text/plain";
				h.add("Content-Type", type);
			}
			
			if(verbose2){System.out.println("Sending: "+response);}
			
			t.sendResponseHeaders(200, response.length());
			OutputStream os = t.getResponseBody();
			os.write(response.getBytes());
			os.close();
		}
		
		public final boolean skipNonCanonical;
	}
	

	
	String toResponse(boolean skipNonCanonical, String[] params){
//		System.err.println("a");
		if(params.length<3){
			if(params.length==2 && "advice".equalsIgnoreCase(params[1])){return TAX_ADVICE;}
			return USAGE;
		}
//		System.err.println("b");
		
		if(params.length>4){return USAGE;}
//		System.err.println("c");
		
		final String query=params[params.length-1];
		final String[] names=query.split(",");
//		System.err.println("d");
		final boolean ancestor=(params.length>3 && params[2].equalsIgnoreCase("ancestor"));
//		System.err.println("e");
//		System.err.println(params[2]+", "+ancestor);
//		System.err.println("f");
		final int type;
		{
			String typeS=params[1];
			Integer value=typeMap.get(typeS);
			if(value==null){
				if(typeS.equalsIgnoreCase("advice")){
					return TAX_ADVICE;
				}else{
					return "{\"error\": \"Bad type; should be gi, taxid, or name.\"}";
				}
			}
			type=value.intValue();
		}
		final int type2=type&15;
		final boolean plaintext=(type>=PT_OFFSET);
		final boolean semicolon=(type>=SC_OFFSET);
		
		if(verbose2){System.out.println("Type: "+type);}
		if(type2==NAME || type2==HEADER){
			for(int i=0; i<names.length; i++){
				if(names[i].contains("%20")){names[i]=names[i].replace("%20", " ");}
				if(type2==HEADER){
					if(names[i].contains("%7C")){names[i]=names[i].replace("%7C", "|");}
					if(names[i].startsWith("@") || names[i].startsWith(">")){names[i]=names[i].substring(1);}
					if(names[i].startsWith("%3E")){names[i]=names[i].substring(3);}
				}
			}
			if(verbose2){System.out.println("Revised: "+Arrays.toString(names));}
		}
		
		if(ancestor){
			if(verbose2){System.out.println("toAncestor: "+Arrays.toString(names));}
			return toAncestor(type, names, plaintext, semicolon, query, skipNonCanonical, !skipNonCanonical);
		}
		
		if(semicolon){
			return toSemicolon(type, names, skipNonCanonical);
		}
		
		if(plaintext){
			return toText(type, names);
		}
		
		if(names.length==1){
			return toJson(type, names[0], skipNonCanonical, !skipNonCanonical).toString();
		}
		ArrayList<JsonObject> list=new ArrayList<JsonObject>();
		for(String name : names){
			list.add(toJson(type, name, skipNonCanonical, !skipNonCanonical));
		}
		return JsonObject.toString(list);
	}
	
	String toAncestor(final int type, final String[] names, boolean plaintext, boolean semicolon, String query, final boolean skipNonCanonical, boolean originalLevel){
		IntList ilist=toIntList(type, names);
		int id=FindAncestor.findAncestor(tree, ilist);
		TaxNode tn=(id>-1 ? tree.getNode(id) : null);
		if(tn==null){return new JsonObject(query, "error","Not found.").toString();}
		if(semicolon){
			return toSemicolon(tn, skipNonCanonical);
		}
		if(plaintext){return ""+id;}
		
		JsonObject j=new JsonObject(query);
		j.add("name", tn.name);
		j.add("tax_id", ""+tn.id);
		j.add("level", ""+tn.levelStringExtended(originalLevel));
		while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=131567){
			if(!skipNonCanonical || tn.isSimple()){
				j.add(toJson(tn, originalLevel));
			}
			if(tn.pid==tn.id){break;}
			tn=tree.getNode(tn.pid);
		}
		return j.toString();
	}
	
	IntList toIntList(final int type, final String[] names){
		IntList list=new IntList(names.length);
		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type2==NAME){
			for(String name : names){
				TaxNode tn=getTaxNodeByName(name);
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type2==NCBI){
			for(String name : names){
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn!=null){list.add(tn.id);}
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				int ncbi=AccessionToTaxid.get(name);
				if(ncbi>=0){list.add(ncbi);}
			}
		}else{
			throw new RuntimeException("{\"error\": \"Bad type\"}");
		}
		return list;
	}
	
	String toText(final int type, final String[] names){
		
		StringBuilder sb=new StringBuilder();
		String comma="";

		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==NCBI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				sb.append(comma);
				int ncbi=AccessionToTaxid.get(name);
				sb.append(ncbi);
				comma=",";
			}
		}else if(type2==HEADER){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeHeader(name);
				if(tn==null){sb.append("-1");}
				else{sb.append(tn.id);}
				comma=",";
			}
		}else{
			return "Bad type; should be pt_gi or pt_name; e.g. /tax/pt_gi/1234";
		}
		
		return sb.toString();
	}
	
	String toSemicolon(final int type, final String[] names, boolean skipNonCanonical){
		
		StringBuilder sb=new StringBuilder();
		String comma="";
		
		int type2=type&15;
		if(type2==GI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeGi(Integer.parseInt(name));
				if(tn==null){sb.append("Not found");}
				else{sb.append(toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==NAME){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeByName(name);
				if(tn==null){sb.append("Not found");}
				else{sb.append(toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==NCBI){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeNcbi(Integer.parseInt(name));
				if(tn==null){sb.append("Not found");}
				else{sb.append(toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==ACCESSION){
			for(String name : names){
				sb.append(comma);
				int ncbi=AccessionToTaxid.get(name);
				TaxNode tn=tree.getNode(ncbi, true);
				if(tn==null){sb.append("Not found");}
				else{sb.append(toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else if(type2==HEADER){
			for(String name : names){
				sb.append(comma);
				TaxNode tn=getTaxNodeHeader(name);
				if(tn==null){sb.append("Not found");}
				else{sb.append(toSemicolon(tn, skipNonCanonical));}
				comma=",";
			}
		}else{
			return "Bad type; should be sc_gi or sc_name; e.g. /tax/sc_gi/1234";
		}
		
		return sb.toString();
	}
	
	String toSemicolon(final TaxNode tn0, boolean skipNonCanonical){
		StringBuilder sb=new StringBuilder();
		if(tn0==null){return "Not found";}
		String semi="";
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(20);
		TaxNode tn=tn0;
		while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=131567){
			if(!skipNonCanonical || tn.isSimple()){
				list.add(tn);
			}
			if(tn.pid==tn.id){break;}
			tn=tree.getNode(tn.pid);
		}
		if(list.isEmpty()){list.add(tn0);}
		boolean addTaxLevel=true;
		for(int i=list.size()-1; i>=0; i--){
			sb.append(semi);
			tn=list.get(i);
			if(addTaxLevel && tn.canonical() && !tn.levelChanged() && tn.isSimple()){
				sb.append(tn.levelToStringShort()).append(':');
			}
			sb.append(tn.name);
			semi=";";
		}
		return sb.toString();
	}
	
	JsonObject toJson(final int type, final String name, boolean skipNonCanonical, boolean originalLevel){
		TaxNode tn=null;
		
		if(type==GI){
			tn=getTaxNodeGi(Integer.parseInt(name));
		}else if(type==NAME){
			tn=getTaxNodeByName(name);
		}else if(type==NCBI){
			tn=getTaxNodeNcbi(Integer.parseInt(name));
		}else if(type==ACCESSION){
			int ncbi=AccessionToTaxid.get(name);
			tn=(ncbi>=0 ? tree.getNode(ncbi) : null);
		}else if(type==HEADER){
			tn=getTaxNodeHeader(name);
		}else{
			return new JsonObject(""+type,"error","Bad type; should be gi, taxid, or name; e.g. /tax/name/homo_sapiens");
		}
		if(verbose2){System.out.println("Got node: "+tn);}
		
		if(tn!=null){
			JsonObject j=new JsonObject(name);
			j.add("name", tn.name);
			j.add("tax_id", ""+tn.id);
			j.add("level", ""+tn.levelStringExtended(originalLevel));
			while(tn!=null && tn.levelExtended!=TaxTree.LIFE_E && tn.id!=131567){
//				System.err.println(tn+", "+(!skipNonCanonical)+", "+tn.isSimple());
				if(!skipNonCanonical || tn.isSimple()){
					j.add(toJson(tn, originalLevel));
//					System.err.println(j);
				}
				if(tn.pid==tn.id){break;}
				tn=tree.getNode(tn.pid);
			}
			return j;
		}
		return new JsonObject(name, "error","Not found.");
	}
	
	JsonObject toJson(TaxNode tn, boolean originalLevel){
		JsonObject j=new JsonObject(tn.levelStringExtended(originalLevel));
		j.add("name",tn.name);
		j.add("tax_id",""+tn.id);
		return j;
	}
	
	TaxNode getTaxNodeByName(String name){
		if(verbose2){System.out.println("Fetching node for "+name);}
		List<TaxNode> list=tree.getNodesByNameExtended(name);
		if(verbose2){System.out.println("Fetched "+list);}
		if(list==null){
			if(verbose2){System.out.println("Fetched in common map "+name);}
			String name2=commonMap.get(name);
			if(verbose2){System.out.println("Fetched "+name2);}
			if(name2!=null){list=tree.getNodesByName(name2);}
		}
		return list==null ? null : list.get(0);
	}
	
	TaxNode getTaxNodeGi(int gi){
		int ncbi=-1;
		try {
			ncbi=GiToNcbi.getID(gi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return ncbi<0 ? null : getTaxNodeNcbi(ncbi);
	}
	
	TaxNode getTaxNodeHeader(String header){
		return tree.getNode(header, '|');
	}
	
	TaxNode getTaxNodeNcbi(int ncbi){
		TaxNode tn=null;
		try {
			tn=tree.getNode(ncbi);
		} catch (Throwable e) {
			if(verbose){e.printStackTrace();}
		}
		return tn;
	}
	
	/** This is called if the program runs with no parameters */
	private void printOptions(){
		throw new RuntimeException("TODO");
	}

	private HashMap<String, Integer> makeTypeMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(63);
		map.put("gi", GI);
		map.put("name", NAME);
		map.put("tax_id", NCBI);
		map.put("ncbi", NCBI);
		map.put("taxid", NCBI);
		map.put("id", NCBI);
		map.put("header", HEADER);
		map.put("accession", ACCESSION);
		map.put("pt_gi", PT_GI);
		map.put("pt_name", PT_NAME);
		map.put("pt_tax_id", PT_NCBI);
		map.put("pt_id", PT_NCBI);
		map.put("pt_ncbi", PT_NCBI);
		map.put("pt_taxid", PT_NCBI);
		map.put("pt_header", PT_HEADER);
		map.put("pt_header", PT_HEADER);
		map.put("pt_accession", PT_ACCESSION);
		map.put("sc_gi", SC_GI);
		map.put("sc_name", SC_NAME);
		map.put("sc_tax_id", SC_NCBI);
		map.put("sc_id", SC_NCBI);
		map.put("sc_ncbi", SC_NCBI);
		map.put("sc_taxid", SC_NCBI);
		map.put("sc_header", SC_HEADER);
		map.put("sc_header", SC_HEADER);
		map.put("sc_accession", SC_ACCESSION);
		
		return map;
	}
	
	public static HashMap<String, String> makeCommonMap(){
		HashMap<String, String> map=new HashMap<String, String>();
		map.put("human", "homo sapiens");
		map.put("cat", "felis catus");
		map.put("dog", "canis lupus familiaris");
		map.put("mouse", "mus musculus");
		map.put("cow", "bos taurus");
		map.put("bull", "bos taurus");
		map.put("horse", "Equus ferus");
		map.put("pig", "Sus scrofa domesticus");
		map.put("sheep", "Ovis aries");
		map.put("goat", "Capra aegagrus");
		map.put("turkey", "Meleagris gallopavo");
		map.put("fox", "Vulpes vulpes");
		map.put("chicken", "Gallus gallus domesticus");
		map.put("wolf", "canis lupus");
		map.put("fruitfly", "drosophila melanogaster");
		map.put("zebrafish", "Danio rerio");
		map.put("catfish", "Ictalurus punctatus");
		map.put("trout", "Oncorhynchus mykiss");
		map.put("salmon", "Salmo salar");
		map.put("tilapia", "Oreochromis niloticus");
		map.put("e coli", "Escherichia coli");
		map.put("e.coli", "Escherichia coli");

		map.put("lion", "Panthera leo");
		map.put("tiger", "Panthera tigris");
		map.put("bear", "Ursus arctos");
		map.put("deer", "Odocoileus virginianus");
		map.put("coyote", "Canis latrans");

		map.put("corn", "Zea mays subsp. mays");
		map.put("maize", "Zea mays subsp. mays");
		map.put("oat", "Avena sativa");
		map.put("wheat", "Triticum aestivum");
		map.put("rice", "Oryza sativa");
		map.put("potato", "Solanum tuberosum");
		map.put("barley", "Hordeum vulgare");
		map.put("poplar", "Populus alba");
		map.put("lettuce", "Lactuca sativa");
		map.put("beet", "Beta vulgaris");
		map.put("strawberry", "Fragaria x ananassa");
		map.put("orange", "Citrus sinensis");
		map.put("lemon", "Citrus limon");
		map.put("soy", "Glycine max");
		map.put("soybean", "Glycine max");
		map.put("grape", "Vitis vinifera");
		map.put("olive", "Olea europaea");
		map.put("cotton", "Gossypium hirsutum");
		map.put("apple", "Malus pumila");
		map.put("bannana", "Musa acuminata");
		map.put("tomato", "Solanum lycopersicum");
		map.put("sugarcane", "Saccharum officinarum");
		map.put("bean", "Phaseolus vulgaris");
		map.put("onion", "Allium cepa");
		map.put("garlic", "Allium sativum");
		
		map.put("pichu", "mus musculus");
		map.put("pikachu", "mus musculus");
		map.put("vulpix", "Vulpes vulpes");
		map.put("ninetails", "Vulpes vulpes");
		map.put("mareep", "Ovis aries");
		
		return map;
	}
	
	/*--------------------------------------------------------------*/
	
	
	
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private String tableFile=TaxTree.defaultTableFile();
	private String treeFile=TaxTree.defaultTreeFile();
	private String accessionFile=null;
	
	private final TaxTree tree;
	private final HashMap<String, Integer> typeMap;
	private final HashMap<String, String> commonMap;
	
	/** Reverse order for tax lines */
	private boolean reverseOrder=true;
	
	/*--------------------------------------------------------------*/
	/*----------------         Final Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	public final int port;
	public final String killCode;
	public String domain="taxonomy.jgi-psf.org";
	
	public final HttpServer httpServer;

	public static final int PT_OFFSET=16;
	public static final int SC_OFFSET=32;
	public static final int UNKNOWN=0, GI=1, NAME=2, NCBI=3, HEADER=4, ACCESSION=5;
	public static final int PT_GI=GI+PT_OFFSET, PT_NAME=NAME+PT_OFFSET, PT_NCBI=NCBI+PT_OFFSET, PT_HEADER=HEADER+PT_OFFSET, PT_ACCESSION=ACCESSION+PT_OFFSET;
	public static final int SC_GI=GI+SC_OFFSET, SC_NAME=NAME+SC_OFFSET, SC_NCBI=NCBI+SC_OFFSET, SC_HEADER=HEADER+SC_OFFSET, SC_ACCESSION=ACCESSION+SC_OFFSET;
	
	public static final String TAX_ADVICE="This site does not give tax advice.";
	public static final String BAD_CODE="Incorrect code.";
	public final String USAGE;
	
	public final SketchSearcher searcher=new SketchSearcher();
	
	
	/*--------------------------------------------------------------*/
	/*----------------        Common Fields         ----------------*/
	/*--------------------------------------------------------------*/
	
	/** Print status messages to this output stream */
	private PrintStream outstream=System.err;
	/** Print verbose messages */
	public static boolean verbose=false, verbose2=false;
	/** True if an error was encountered */
	public boolean errorState=false;
	
}
