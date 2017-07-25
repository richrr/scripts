package tax;

import java.io.PrintStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.regex.Pattern;

import fileIO.ReadWrite;
import fileIO.TextFile;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.IntList;

/**
 * @author Brian Bushnell
 * @date Mar 6, 2015
 *
 */
public class TaxTree implements Serializable{
	
	private static final long serialVersionUID = -5423000353519543015L;
	
	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(String[] args){
		ReadWrite.USE_UNPIGZ=true;
		ReadWrite.USE_PIGZ=true;
		ReadWrite.ZIPLEVEL=(Shared.threads()>2 ? 11 : 9);
		ReadWrite.PIGZ_BLOCKSIZE=256;
		ReadWrite.PIGZ_ITERATIONS=60;
		for(String arg : args){
			String[] split=arg.split("=");
			String a=split[0].toLowerCase();
			String b=split.length>1 ? split[1] : null;
			if(b==null || b.equalsIgnoreCase("null")){b=null;}
			dna.Parser.parseZip(arg, a, b);
		}
		Timer t=new Timer();
		TaxTree tree=new TaxTree(args[0], args[1]);
		t.stop();
		System.out.println("Retained "+tree.nodeCount+" nodes:");
		for(int i=tree.treeLevelsExtended.length-1; i>=0; i--){
			System.out.print(tree.nodesPerLevelExtended[i]+"\t"+taxaNamesExtended[i]);
			if(verbose){
				int lim=10;
				for(int j=0; j<lim && j<tree.treeLevelsExtended[i].length; j++){
					TaxNode n=tree.treeLevelsExtended[i][j];
					System.out.print("\n"+n+" -> "+tree.nodes[n.pid]);
				}
				for(int j=tree.treeLevelsExtended[i].length-lim; j<tree.treeLevelsExtended[i].length; j++){
					if(j>=lim){
						TaxNode n=tree.treeLevelsExtended[i][j];
						System.out.print("\n"+n+" -> "+tree.nodes[n.pid]);
					}
				}
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("Time: \t"+t);
		
		if(args.length>2){//Write a tree
			ReadWrite.write(tree, args[2], true);
		}
	}
	
	public TaxTree(String namesFile, String nodesFile){
		
		nodes=getNames(namesFile);
		getNodes(nodesFile, nodes);
		
		if(simplify){
			if(verbose){System.out.println("Simplifying.");}
			int removed=simplify(nodes);
			if(verbose){System.out.println("Removed "+removed+" nodes.");}
		}
		assert(test(nodes)==0);
		
		for(TaxNode n : nodes){
			if(n!=null){
				nodesPerLevel[n.level]++;
				nodesPerLevelExtended[n.levelExtended]++;
			}
		}
		
//		for(int i=0; i<nodesPerLevel.length; i++){
//			treeLevels[i]=new TaxNode[nodesPerLevel[i]];
//		}
		for(int i=0; i<nodesPerLevelExtended.length; i++){
			treeLevelsExtended[i]=new TaxNode[nodesPerLevelExtended[i]];
		}
		
//		{
//			int[] temp=new int[nodesPerLevel.length];
//			for(TaxNode n : nodes){
//				if(n!=null){
//					int level=n.level;
//					treeLevels[level][temp[level]]=n;
//					temp[level]++;
//				}
//			}
//		}
		
		{
			int[] temp=new int[nodesPerLevelExtended.length];
			for(TaxNode n : nodes){
				if(n!=null){
					int level=n.levelExtended;
					treeLevelsExtended[level][temp[level]]=n;
					temp[level]++;
				}
			}
		}
		nodeCount=(int)Tools.sum(nodesPerLevelExtended);
		
	}
	
	public static final TaxTree loadTaxTree(String taxTreeFile, PrintStream outstream, boolean hashNames){
		if(taxTreeFile==null){return null;}
		return loadTaxTree(taxTreeFile, null, null, outstream, hashNames);
	}
	
	public static final TaxTree loadTaxTree(String taxTreeFile, String taxNameFile, String taxNodeFile, PrintStream outstream, boolean hashNames){
		assert(taxTreeFile!=null || (taxNameFile!=null && taxNodeFile!=null)) : "Must specify both taxname and taxnode files.";
		Timer t=new Timer();
		if(outstream!=null){outstream.print("\nLoading tax tree; ");}
		final TaxTree tree;
		if(taxTreeFile!=null){
			tree=ReadWrite.read(TaxTree.class, taxTreeFile, true);
		}else{
			tree=new TaxTree(taxNameFile, taxNodeFile);
		}
		t.stop();
		if(hashNames){
			outstream.println("Hashing names.");
			tree.hashNames();
		}
		if(outstream!=null){
			outstream.println("time: \t"+t);
			Shared.printMemory();
			outstream.println();
		}
		return tree;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------         Construction         ----------------*/
	/*--------------------------------------------------------------*/
	
	private static TaxNode[] getNames(String fname){
		ArrayList<TaxNode> list=new ArrayList<TaxNode>(200000);
		int max=0;
		
		TextFile tf=new TextFile(fname, false, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			if(s.contains("scientific name")){
				String[] split=delimiter.split(s, 3);
				assert(split.length==3) : s;
				int id=Integer.parseInt(split[0]);
				String name=split[1];
				if(id==1 && name.equalsIgnoreCase("root")){name="Life";}
				max=Tools.max(max, id);
				list.add(new TaxNode(id, name));
			}
		}
		
		TaxNode[] nodes=new TaxNode[max+1];
		for(TaxNode n : list){
			assert(nodes[n.id]==null || nodes[n.id].equals(n)) : nodes[n.id]+" -> "+n;
			nodes[n.id]=n;
		}
		
		return nodes;
	}
	
	public void hashNames(){
		assert(nameMap==null);
		assert(nameMapLower==null);
		nameMap=new HashMap<String, ArrayList<TaxNode>>((int)Tools.mid(2, nodes.length*1.5, Integer.MAX_VALUE));
		nameMapLower=new HashMap<String, ArrayList<TaxNode>>((int)Tools.mid(2, nodes.length*1.5, Integer.MAX_VALUE));
		for(TaxNode n : nodes){
			if(n!=null){
				String name=n.name;
				if(name.indexOf('_')>=0){
					name=name.replace('_', ' ').trim();
				}
				if(name!=null && !name.equals("environmental samples")){
					{
						ArrayList<TaxNode> list=nameMap.get(name);
						if(list==null){
							list=new ArrayList<TaxNode>();
							nameMap.put(name, list);
						}
						list.add(n);
					}
					{
						String lc=name.toLowerCase();
						ArrayList<TaxNode> list=nameMapLower.get(lc);
						if(list==null){
							list=new ArrayList<TaxNode>();
							nameMapLower.put(lc, list);
						}
						list.add(n);
					}
				}
			}
		}
	}
	
	private static TaxNode[] getNodes(String fname, TaxNode[] nodes){
		
		int max=0;
		
		LinkedHashMap<String, int[]> oddNames=new LinkedHashMap<String, int[]>();
		
		TextFile tf=new TextFile(fname, false, false);
		for(String s=tf.nextLine(); s!=null; s=tf.nextLine()){
			String[] split=delimiter.split(s, 4);
			assert(split.length==4) : s;
			int id=-1, pid=-1, level=-1, levelExtended=-1;
			
			id=Integer.parseInt(split[0]);
			try {
				pid=Integer.parseInt(split[1]);
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.err.println("Bad line: "+s+"\n"+Arrays.toString(split));
			}
			boolean alt=false;
			{
				String key=split[2];
				Integer obj0=levelMap.get(key);
				Integer obj=levelMapExtended.get(key);
				if(obj0==null){
					obj0=altLevelMap.get(key);
					alt=true;
				}
				if(obj0!=null){
					level=obj0;
					levelExtended=obj;
					if(id==pid){
						level=LIFE;
						levelExtended=LIFE_E;
						alt=false;
					}
				}else{
					if(id==pid){
						level=LIFE;
						levelExtended=LIFE_E;
						alt=false;
					}else{
						int[] count=oddNames.get(key);
						if(count==null){
							count=new int[1];
							oddNames.put(key, count);
						}
						count[0]++;
					}
				}
			}
			max=Tools.max(max, id);
			TaxNode n=nodes[id];
			assert(n!=null && n.pid<0) : n+" -> "+s;
			n.pid=pid;
			n.level=level;
			n.levelExtended=levelExtended;
			n.setOriginalLevel(levelExtended);
			n.setCanonical(!alt);
			assert(n.canonical()==n.isSimple() || n.levelExtended==NO_RANK_E) : n.canonical()+", "+n.isSimple()+", "+n.level+", "+n.levelExtended;
		}
		
		if(oddNames.size()>0){
			System.out.println("Found "+oddNames.size()+" unknown taxonomic levels:");
			if(verbose){
				for(String s : oddNames.keySet()){
					System.out.println(oddNames.get(s)[0]+"\t"+s);
				}
			}
		}
		
		return nodes;
	}
	
	private int simplify(TaxNode nodes[]){
		
		int failed=test(nodes);
		
		int removed=0;
		int reassigned=0;
		
		if(reassign){
			boolean changed=true;
			int changedCount=0;
			while(changed){
				changed=false;
				for(int i=0; i<nodes.length; i++){
					TaxNode n=nodes[i];
					if(nodes[i]!=null && nodes[i].levelExtended<1){
						int pid=n.pid;
						TaxNode parent=nodes[pid];
						assert(parent!=null) : n;
						if(parent.levelExtended==SPECIES_E || parent.levelExtended==SUBSPECIES_E){
//							assert(n.id!=594){
//								System.err.println(n+"->"+parent+"\n")
//							}
							changed=true;
							n.levelExtended=Tools.max(SUBSPECIES_E, parent.levelExtended-1);
							n.level=Tools.max(SUBSPECIES, parent.level-1);
							changedCount++;
						}
					}
				}
			}
			System.err.println("Assigned levels to "+changedCount+" unranked nodes.");
		}
		
		
		if(skipNorank){//Skip nodes with unknown taxa
			if(verbose){System.out.println("A0");}
			
			for(int i=0; i<nodes.length; i++){
				TaxNode n=nodes[i];
				if(n!=null){
					int pid=n.pid;
					TaxNode parent=nodes[pid];
					assert(parent!=null) : n;
					assert(parent!=n || pid==1) : n+", "+pid;
					while(parent.levelExtended<1 && n.levelExtended>parent.levelExtended){
						//System.err.println("Reassigned from "+parent);
						assert(parent.id!=parent.pid);
						parent=nodes[parent.pid];
						n.pid=parent.id;
						reassigned++;
					}
				}
			}
			
			for(int i=0; i<nodes.length; i++){
				if(nodes[i]!=null && nodes[i].levelExtended<0){
					System.err.println("Removed "+nodes[i]);
					nodes[i]=null;
					removed++;
				}
			}
			if(verbose){System.out.println("Skipped "+reassigned+" unranked parents, removed "+removed+" invalid nodes.");}
		}
		
		if(inferRankLimit>0){//Infer level for unset nodes (from "no rank")
			if(verbose){System.out.println("A");}
			int changed=1;
			while(changed>0){
				changed=0;
				for(final TaxNode n : nodes){
					if(n!=null){
						if(n.levelExtended==0){
							TaxNode parent=nodes[n.pid];
							if(n!=parent && parent.levelExtended>0 && parent.levelExtended<=inferRankLimit+1){
								n.levelExtended=Tools.max(1, parent.levelExtended-1);
								assert(n.levelExtended>0 && n.levelExtended<=parent.levelExtended && n.levelExtended<=inferRankLimit);
								changed++;
							}
						}
					}
				}
				if(verbose){System.out.println("changed: "+changed);}
			}
			
//			System.out.println("B");
//			for(TaxNode n : nodes){
//				if(n!=null && n.level==0){
//					n.level=-1;
//				}
//			}
		}
		
		failed=test(nodes);
		
//		if(reassign){//Skip nodes with duplicate taxa
//			if(verbose){System.out.println("D");}
//			int changed=1;
//			while(changed>0){
//				changed=0;
//				for(final TaxNode n : nodes){
//					if(n!=null){
//						TaxNode parent=nodes[n.pid];
//						TaxNode grandparent=nodes[parent.pid];
//						assert(n.level<=parent.level || parent.level<1 || !parent.canonical()) : n+" -> "+parent+" -> "+grandparent;
//						assert(parent.level<=grandparent.level || grandparent.level<1 || !grandparent.canonical()) : n+" -> "+parent+" -> "+grandparent;
//
//						while(parent!=grandparent && (parent.level<0 || (parent.level==grandparent.level && !parent.canonical()) || 
//								n.level>parent.level || (n.level==parent.level))){
//							parent=grandparent;
//							grandparent=nodes[parent.pid];
//							n.pid=parent.id;
//							reassigned++;
//							changed++;
//						}
//					}
//				}
//				if(verbose){System.out.println("changed: "+changed);}
//			}
//			if(verbose){System.out.println("E");}
//			for(int i=0; i<nodes.length; i++){
//				if(nodes[i]!=null && nodes[i].level<0){
//					nodes[i]=null;
//					removed++;
//				}
//			}
//		}
		
		failed=test(nodes);

		if(verbose){System.out.println("F");}
		{//Ensure the tree is now clean
			for(int i=0; i<nodes.length; i++){
				TaxNode n=nodes[i];
				if(n!=null){
					TaxNode parent=nodes[n.pid];
					TaxNode grandparent=nodes[parent.pid];
					assert(n==parent || n.levelExtended<=parent.levelExtended || !n.canonical() || n.levelExtended<1 || parent.levelExtended<1) : n+" -> "+parent+" -> "+grandparent;
					assert(parent==grandparent || parent.levelExtended<=grandparent.levelExtended || !parent.canonical() || parent.levelExtended<1 || grandparent.levelExtended<1) : n+" -> "+parent+" -> "+grandparent;
				}
			}
		}
		
//		if(verbose){System.err.println("Reassignments: "+reassigned);}
		
		return removed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------          Validation          ----------------*/
	/*--------------------------------------------------------------*/
	
	private static int test(TaxNode[] nodes){
		int failed=0;
		for(final TaxNode n : nodes){
			if(n!=null){
				TaxNode parent=nodes[n.pid];
				assert(n==parent || n.level<=parent.level || parent.level<1 || !parent.canonical()) : n+" -> "+parent;
				assert(n==parent || n.levelExtended<=parent.levelExtended || parent.levelExtended<1) : n+" -> "+parent;
//				assert(n==parent || n.level<parent.level || parent.level<1 || !n.canonical() || !parent.canonical()) : n+" -> "+parent;
				if(n!=parent && n.level>parent.level && parent.level>=1 && n.canonical() && parent.canonical()){
					if(verbose){System.out.println("Error: "+n+" -> "+parent);}
					failed++;
				}else if(n!=parent && parent.levelExtended>=1 && n.levelExtended>=parent.levelExtended){
//					if(verbose){System.out.println("Error: "+n+" -> "+parent);}
//					failed++;
				}
				assert(n!=parent || n.id<=1) : n;
			}
		}
		if(verbose){System.out.println(failed+" nodes failed.");}
		return failed;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------        Outer Methods         ----------------*/
	/*--------------------------------------------------------------*/
	
	
	public static int getID(String s){return GiToNcbi.getID(s);}
	
	public static int getID(byte[] s){return GiToNcbi.getID(s);}
	
	/** Return the ancestor with taxonomic level at least minLevel */
	public TaxNode getNode(String s, int minLevel){
		TaxNode tn=getNode(s, '|');
		final int minLevelExtended=levelToExtended(minLevel);
		while(tn!=null && tn.levelExtended<minLevelExtended && tn.pid!=tn.id){
			tn=getNode(tn.pid);
		}
		return tn;
	}
	
	public int commonAncestor(final int a, final int b){
		TaxNode an=getNode(a), bn=getNode(b);
		assert(an!=null) : "Invalid taxID: "+a;
		assert(bn!=null) : "Invalid taxID: "+b;
		TaxNode cn=commonAncestor(an, bn);
		assert(cn!=null) : "No common ancestor: "+an+", "+bn;
		if(cn==null){return -1;}
		return cn.id;
	}
	
	public TaxNode commonAncestor(TaxNode a, TaxNode b){
		assert(a!=null && b!=null) : "Null parameters.";
		if(a==null){return b;}
		if(b==null){return a;}
		if(a==null && b==null){return null;}
		
		while(a!=b){
			if(a.levelExtended<b.levelExtended){
				a=getNode(a.pid);
			}else{
				b=getNode(b.pid);
			}
		}
		return a;
	}
	
	public TaxNode getNode(String s, char delimiter){
		{
			int index=s.indexOf(delimiter);
			if(index<0){
				delimiter='~';
				index=s.indexOf(delimiter);
				if(index<0){
					delimiter='_';
					index=s.indexOf(delimiter);
				}
			}
			int number=-1;
			
			Throwable e=null;
			
			if(index==2 && s.length()>3 && s.startsWith("gi") && Character.isDigit(s.charAt(3))){
//				System.err.println("Parsing gi number.");
				try {
					number=GiToNcbi.parseGiToNcbi(s, delimiter);
				} catch (Throwable e2) {
					e=e2;
				}
//				if(number!=-1){System.err.println("number="+number);}
			}else if(index==4 && s.length()>5 && s.startsWith("ncbi") && Character.isDigit(s.charAt(5))){
//				System.err.println("Parsing ncbi number.");
				number=GiToNcbi.getID(s);
			}
			
			if(number<0 && index>=0 && (delimiter=='|' || delimiter=='~')){
				String[] split=(delimiter=='|' ? delimiterPipe.split(s) : delimiterTilde.split(s));
				if(AccessionToTaxid.LOADED){
					number=parseAccessionToTaxid(split);	
				}
				if(number<0){
					number=parseHeaderNameToTaxid(split);
				}
			}
			
			if(number<0 && e!=null){
				assert(false) : e;
				throw new RuntimeException(e);
			}
			
			//TaxServer code could go here...
			
			if(number>=0){return getNode(number);}
		}
		if(verbose){System.err.println("Can't process name "+s);}
		if(Character.isDigit(s.charAt(0)) && s.length()<=9){
			try {
				return getNode(Integer.parseInt(s));
			} catch (NumberFormatException e) {
				//ignore
			}
		}
		return null;
	}
	
	public int parseAccessionToTaxid(String[] split){
		if(split.length<4){
			return -1;
		}
		int ncbi=AccessionToTaxid.get(split[3]);
		return ncbi;
	}
	
	public int parseHeaderNameToTaxid(String[] split){
		if(split.length<5){
			return -1;
		}
		return parseNameToTaxid(split[4]);
	}
	
	public int parseNameToTaxid(String name){
		List<TaxNode> list=null;
		
		list=getNodesByNameExtended(name);
		
		if(list==null || list.size()>1){return -1;}
		return list.get(0).id;
	}
	
	public List<TaxNode> getNodesByNameExtended(String name){
		List<TaxNode> list=null;
		
		list=getNodesByName(name);
		if(list!=null){return list;}
		
		name=name.replaceAll("_", " ").trim();
		list=getNodesByName(name);
		if(list!=null){return list;}
		
		String[] split2=name.split(" ");
		
		if(split2.length>7){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6]+" "+split2[7];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		if(split2.length>6){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5]+" "+split2[6];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		if(split2.length>5){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4]+" "+split2[5];
			list=getNodesByName(term);
//			System.err.println("6:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>4){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3]+" "+split2[4];
			list=getNodesByName(term);
//			System.err.println("5:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>3){
			String term=split2[0]+" "+split2[1]+" "+split2[2]+" "+split2[3];
			list=getNodesByName(term);
//			System.err.println("4:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>2){
			String term=split2[0]+" "+split2[1]+" "+split2[2];
			list=getNodesByName(term);
//			System.err.println("3:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>1){
			String term=split2[0]+" "+split2[1];
			list=getNodesByName(term);
//			System.err.println("2:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		if(split2.length>0){
			String term=split2[0];
			list=getNodesByName(term);
//			System.err.println("1:\n"+Arrays.toString(split)+"\n"+Arrays.toString(split2)+"\n"+term+" -> "+list);
			if(list!=null){return list;}
		}
		
		return null;
	}
	
//	public TaxNode getNode(byte[] s){
//		if(Tools.indexOf(s, (byte)'|')>=0){return getNode(GiToNcbi.getID(s));}
//		
//		{
//			int index=Tools.indexOf(s, (byte)'|');
//			if(index<0){index=Tools.indexOf(s, (byte)'_');}
//			int number=-1;
//			if(index==2 && s.length>3 && Tools.startsWith(s, "gi") && Character.isDigit(s[3])){
////				System.err.println("Parsing gi number.");
//				number=GiToNcbi.parseGiToNcbi(s);
//			}else if(index==4 && s.length>5 && Tools.startsWith(s, "ncbi") && Character.isDigit(s[5])){
////				System.err.println("Parsing ncbi number.");
//				number=GiToNcbi.getID(s);
//			}
//			if(number>=0){return getNode(number);}
//		}
//		if(verbose){System.err.println("Can't process name "+new String(s));}
//		if(Character.isDigit(s[0]) && s.length<=9){
//			try {
//				return getNode(Tools.parseInt(s, 0, s.length));
//			} catch (NumberFormatException e) {
//				//ignore
//			}
//		}
//		return null;
//	}
	
	public TaxNode getNode(int id){
		assert(id<nodes.length) : id+", "+nodes.length;
		return id<0 ? null : nodes[id];
	}
	
	public TaxNode getNode(int id, boolean skipAssertion){
		assert(skipAssertion || id<nodes.length) : id+", "+nodes.length;
		return id<0 || id>=nodes.length ? null : nodes[id];
	}
	
	public TaxNode getNodeAtLevel(int id, int minLevel){
		return getNodeAtLevel(id, minLevel, DOMAIN);
	}
	
	public TaxNode getNodeAtLevel(int id, int minLevel, int maxLevel){
		TaxNode tn=getNode(id);
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
		while(tn!=null && tn.pid!=tn.id && tn.levelExtended<minLevelExtended){
			TaxNode temp=getNode(tn.pid);
			if(temp==null || temp.levelExtended>maxLevelExtended){break;}
			tn=temp;
		}
		return tn;
	}

	public TaxNode getNodeByName(String s){
		List<TaxNode> list=getNodesByName(s, false);
		if(list==null){list=getNodesByName(s, true);}
		if(list==null || list.size()<1){return null;}
		if(list.size()==1){return list.get(0);}
		assert(false) : "Found multiple nodes for '"+s+"':\n"+list+"\n";
		return list.get(0);
	}
	public List<TaxNode> getNodesByName(String s){
		List<TaxNode> list=getNodesByName(s, false);
		if(list==null){list=getNodesByName(s, true);}
		return list;
	}
	private List<TaxNode> getNodesByName(String s, boolean lowercase){
		if(s.indexOf('_')>=0){s=s.replace('_', ' ');}
		if(lowercase){s=s.toLowerCase();}
//		System.err.println("Searching for "+s);
		final HashMap<String, ArrayList<TaxNode>> map=(lowercase ? nameMapLower : nameMap);
		ArrayList<TaxNode> list=map.get(s);
		if(list!=null){return list;}
//		System.err.println("No matches for '"+s+"'");
		
//		assert(false) : nameMap.containsKey(s)+", "+nameMapLower.containsKey(s);
		
		if(s.indexOf('_')<0 && s.indexOf(' ')<0){return null;}
		String[] split=delimiter2.split(lowercase ? s.toLowerCase() : s, 8);
//		System.err.println("Array: "+Arrays.toString(split));
		list=map.get(split[split.length-1]);
		if(list==null){return list;}
//		System.err.println(list==null ? "No matches for "+split[split.length-1] : "Found list( "+list.size()+")");
		
		int matchCount=0;
		for(TaxNode tn : list){
			if(tn.matchesName(split, split.length-1, this)){matchCount++;}
		}
		if(matchCount==list.size()){return list;}
		if(matchCount<1){return null;}
		ArrayList<TaxNode> hits=new ArrayList<TaxNode>(matchCount);
		for(TaxNode tn : list){
			if(tn.matchesName(split, split.length-1, this)){hits.add(tn);}
		}
		return hits;
	}
	public ArrayList<TaxNode> getAncestors(int id){
		TaxNode current=getNode(id);
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		while(current!=null && current.pid!=current.id){//ignores root
			list.add(current);
			current=getNode(current.pid);
		}
		//optionally add root here
		return list;
	}
	
	public void increment(IntList ids, IntList counts, boolean sync){
		
		ids.sort();
		ids.getUniqueCounts(counts);
		
		if(!sync){
			for(int i=0; i<ids.size; i++){
				int id=ids.get(i);
				int count=counts.get(i);
				incrementRaw(id, count);
			}
		}else{
			synchronized(this){
				for(int i=0; i<ids.size; i++){
					int id=ids.get(i);
					int count=counts.get(i);
					incrementRaw(id, count);
				}
			}
		}
	}
	
	public void incrementRaw(int id, long amt){
		nodes[id].incrementRaw(amt);
	}
	
	public void percolateUp(){
		for(int i=0; i<treeLevelsExtended.length; i++){
			percolateUp(i);
		}
	}
	
	public void percolateUp(final int fromLevel){
		final TaxNode[] stratum=treeLevelsExtended[fromLevel];
		for(final TaxNode n : stratum){
			n.incrementSum(n.countRaw);
			TaxNode parent=nodes[n.pid];
			if(n!=parent){
				parent.incrementSum(n.countSum);
			}
		}
	}
	
	/** Add this amount to the node and all its ancestors. */
	public void percolateUp(TaxNode node, long amt){
		if(amt==0){return;}
		if(verbose){System.err.println("percolateUp("+amt+") node: "+node);}
		while(node.id!=node.pid){
			node.incrementSum(amt);
			node=nodes[node.pid];
		}
		node.incrementSum(amt);
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit){
		return gatherNodesAtLeastLimit(limit, 0, nodesPerLevelExtended.length-1);
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final long limit, final int minLevel, final int maxLevel){
		final int minLevelExtended=levelToExtended(minLevel);
		final int maxLevelExtended=levelToExtended(maxLevel);
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		for(int i=minLevelExtended; i<nodesPerLevelExtended.length && i<=maxLevelExtended; i++){
			list.addAll(gatherNodesAtLeastLimit(i, limit));
		}
		Shared.sort(list, TaxNode.countComparator);
		return list;
	}
	
	public ArrayList<TaxNode> gatherNodesAtLeastLimit(final int fromLevel, final long limit){
		final int fromLevelExtended=levelToExtended(fromLevel);
		ArrayList<TaxNode> list=new ArrayList<TaxNode>();
		final TaxNode[] stratum=treeLevelsExtended[fromLevelExtended];
		for(final TaxNode n : stratum){
			if(n.countSum>=limit){
				list.add(n);
				TaxNode parent=nodes[n.pid];
				if(n!=parent){
					percolateUp(parent, -n.countSum);
				}
			}
		}
		Shared.sort(list, TaxNode.countComparator);
		return list;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------     Static Initializers      ----------------*/
	/*--------------------------------------------------------------*/

	/**
	 * @return
	 */
	private static HashMap<String, Integer> makeLevelMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(31);
		for(int i=0; i<taxaNames.length; i++){
			map.put(taxaNames[i], i);
			map.put(taxaNames[i].toUpperCase(), i);
		}
		return map;
	}

	/**
	 * @return
	 */
	private static HashMap<String, Integer> makeLevelMapExtended() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(129);
		for(int i=0; i<taxaNamesExtended.length; i++){
			map.put(taxaNamesExtended[i], i);
			map.put(taxaNamesExtended[i].toUpperCase(), i);
		}
		return map;
	}

	/**
	 * @return
	 */
	private static HashMap<String, Integer> makeAltLevelMap() {
		HashMap<String, Integer> map=new HashMap<String, Integer>(67);
		for(int i=0; i<taxaNames.length; i++){
			map.put(taxaNames[i], i);
		}
		
		//Add synonyms
		map.put("subfamily", map.get("family"));
		map.put("tribe", map.get("family"));
		map.put("varietas", map.get("subspecies"));
		map.put("subgenus", map.get("genus"));
		map.put("forma", map.get("subspecies"));
		map.put("species group", map.get("genus"));
		map.put("subclass", map.get("class"));
		map.put("species subgroup", map.get("species"));
		map.put("infraorder", map.get("order"));
		map.put("superorder", map.get("class"));
		map.put("subphylum", map.get("phylum"));
		map.put("infraclass", map.get("class"));
		map.put("superkingdom", map.get("division"));
		map.put("parvorder", map.get("order"));
		map.put("superclass", map.get("phylum"));
		map.put("superphylum", map.get("kingdom"));
		map.put("subkingdom", map.get("kingdom"));
		map.put("superfamily", map.get("order"));
		map.put("superkingdom", map.get("domain"));
		map.put("suborder", map.get("order"));
		map.put("subtribe", map.get("family"));
//		map.put("no rank", map.get("subspecies"));
		
		return map;
	}
	
	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/
	
	public final TaxNode[] nodes;
	public final int[] nodesPerLevel=new int[taxaNames.length];
	public final int[] nodesPerLevelExtended=new int[taxaNamesExtended.length];
	public final int nodeCount;

//	public final TaxNode[][] treeLevels=new TaxNode[taxaNames.length][];
	public final TaxNode[][] treeLevelsExtended=new TaxNode[taxaNamesExtended.length][];
	
	HashMap<String, ArrayList<TaxNode>> nameMap;
	HashMap<String, ArrayList<TaxNode>> nameMapLower;
	
	public int minValidTaxa=0;
	
	public boolean simplify=true;
	public boolean reassign=true;
	public boolean skipNorank=false;
	public int inferRankLimit=0;//levelMap.get("species");
	
	/*--------------------------------------------------------------*/
	/*----------------          Constants           ----------------*/
	/*--------------------------------------------------------------*/

	public static final int stringToLevel(String s){return levelMap.containsKey(s) ? levelMap.get(s) : altLevelMap.get(s);}
	public static final int stringToLevelExtended(String s){return levelMapExtended.get(s);}
	public static final String levelToString(int x){return taxaNames[x];}
	public static final String levelToStringExtended(int x){return taxaNamesExtended[x];}
	public static final String levelToStringShort(int x){return taxaNamesShort[x];}
	
	private static final String[] taxaNames=new String[] {
		"no rank", "subspecies", "species", "genus",
		"family", "order", "class", "phylum",
		"kingdom", "superkingdom", "domain", "life"
	};
	
	private static final String[] taxaNamesExtended=new String[] {
		"no rank", 
		"forma", "varietas", "subspecies",
		"species",
		"species subgroup", "species group", "subgenus", "genus",
		"subtribe", "tribe", "subfamily", "family",
		"superfamily", "parvorder", "infraorder", "suborder", "order",
		"superorder", "infraclass", "subclass", "class", 
		"superclass", "subphylum", "phylum",
		"superphylum", "subkingdom", "kingdom", 
		"superkingdom", "domain", 
		"life"
	};
	
	
	private static final String[] taxaNamesShort=new String[] {
			"nr", "ss", "s", "g",
			"f", "o", "c", "p",
			"k", "sk", "d", "l"
	};
	
	public static final int NO_RANK=0, SUBSPECIES=1, SPECIES=2, GENUS=3,
			FAMILY=4, ORDER=5, CLASS=6, PHYLUM=7, KINGDOM=8, SUPERKINGDOM=9, DOMAIN=10, LIFE=11;
	
	
	private static final HashMap<String, Integer> levelMap=makeLevelMap();
	private static final HashMap<String, Integer> levelMapExtended=makeLevelMapExtended();
	private static final HashMap<String, Integer> altLevelMap=makeAltLevelMap();
	
	public static final int NO_RANK_E=NO_RANK, SUBSPECIES_E=stringToLevelExtended("subspecies"), 
			SPECIES_E=stringToLevelExtended("species"), GENUS_E=stringToLevelExtended("genus"), 
			FAMILY_E=stringToLevelExtended("family"), ORDER_E=stringToLevelExtended("order"),
			CLASS_E=stringToLevelExtended("class"), PHYLUM_E=stringToLevelExtended("phylum"), 
			KINGDOM_E=stringToLevelExtended("kingdom"), SUPERKINGDOM_E=stringToLevelExtended("superkingdom"), 
			DOMAIN_E=stringToLevelExtended("domain"), LIFE_E=stringToLevelExtended("life");

	private static final int[] levelToExtended=new int[] {
			NO_RANK_E, SUBSPECIES_E, SPECIES_E, GENUS_E, FAMILY_E,
			ORDER_E, CLASS_E, PHYLUM_E, KINGDOM_E, SUPERKINGDOM_E, DOMAIN_E, LIFE_E
		};
	
	public static final int levelToExtended(int level){
		return level<0 ? level : levelToExtended[level];
	}
	
	private static final Pattern delimiter = Pattern.compile("\t\\|\t");
	private static final Pattern delimiterPipe = Pattern.compile("\\|");
//	private static final Pattern delimiterBang = Pattern.compile("\\!");
	private static final Pattern delimiterTilde = Pattern.compile("\\~");
	private static final Pattern delimiter2 = Pattern.compile("[\\s_]+");

	public static String TAX_PATH="/global/projectb/sandbox/gaag/bbtools/tax/latest";

	public static final String defaultTableFile(){return defaultTableFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultTreeFile(){return defaultTreeFile.replaceAll("TAX_PATH", TAX_PATH);}
	public static final String defaultAccessionFile(){return defaultAccessionFile.replaceAll("TAX_PATH", TAX_PATH);}
	
	private static final String defaultTableFile="TAX_PATH/gitable.int1d.gz";
	private static final String defaultTreeFile="TAX_PATH/tree.taxtree.gz";
	private static final String defaultAccessionFile="TAX_PATH/dead_nucl.accession2taxid.gz,"
			+ "TAX_PATH/dead_prot.accession2taxid.gz,TAX_PATH/dead_wgs.accession2taxid.gz,"
			+ "TAX_PATH/nucl_est.accession2taxid.gz,TAX_PATH/nucl_gb.accession2taxid.gz,"
			+ "TAX_PATH/nucl_gss.accession2taxid.gz,TAX_PATH/nucl_wgs.accession2taxid.gz,"
			+ "TAX_PATH/pdb.accession2taxid.gz,TAX_PATH/prot.accession2taxid.gz";

	public static boolean verbose=false;
	public static boolean SHOW_WARNINGS=false;
	
}
