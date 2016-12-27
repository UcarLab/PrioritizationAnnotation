package nearestgene;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.sql.Connection;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

public class NearestTSSUtil {
	private TreeMap<String, TSS[]> _tsssorted;

	public NearestTSSUtil(Connection conn, String database, String genecol, String chrcol, String startcol, String endcol, String strandcol) throws SQLException{
		String startsortedsql = "SELECT g."+genecol+", g."+chrcol+", CASE g."+strandcol+" WHEN '-' THEN g."+endcol+" ELSE g."+startcol+" END AS TSS, g."+strandcol+" FROM "+database+" AS g ORDER BY tss ASC";
		PreparedStatement ps = conn.prepareStatement(startsortedsql);
		ResultSet rs = ps.executeQuery();
		
		TreeMap<String, LinkedList<TSS>> map = new TreeMap<String, LinkedList<TSS>>();
		while(rs.next()){
			String chr = rs.getString(2);
			if(!map.containsKey(chr)){
				map.put(chr, new LinkedList<TSS>());
			}
			map.get(chr).add(new TSS(rs.getString(1), chr, rs.getString(4), rs.getInt(3)));
		}
		rs.close();
		ps.close();
		_tsssorted = new TreeMap<String, TSS[]>();
		for(Iterator<Entry<String, LinkedList<TSS>>> it = map.entrySet().iterator(); it.hasNext();){
			Entry<String, LinkedList<TSS>> e = it.next();
			_tsssorted.put(e.getKey(), e.getValue().toArray(new TSS[0]));
		}
		
	}
	
	public NearestTSSUtil(String refflatfile) throws SQLException, IOException{
		TreeMap<String, LinkedList<TSS>> map = parseRefflatFile(refflatfile);
		_tsssorted = new TreeMap<String, TSS[]>();
		for(Iterator<Entry<String, LinkedList<TSS>>> it = map.entrySet().iterator(); it.hasNext();){
			Entry<String, LinkedList<TSS>> e = it.next();
			_tsssorted.put(e.getKey(), sortByTSS(e.getValue()).toArray(new TSS[0]));
		}
	}
	
	public TSS[] getTSSLocations(){
		LinkedList<TSS> rv = new LinkedList<TSS>();
		for(Iterator<Entry<String, TSS[]>> it = _tsssorted.entrySet().iterator(); it.hasNext();){
			Entry<String, TSS[]> next = it.next();
			TSS[] tsslocations = next.getValue();
			for(int i = 0; i < tsslocations.length; i++){
				rv.add(tsslocations[i]);
			}
		}
		return rv.toArray(new TSS[0]);
	}
	
	private TreeMap<String, LinkedList<TSS>> parseRefflatFile(String file) throws IOException{
		TreeMap<String, LinkedList<TSS>> map = new TreeMap<String, LinkedList<TSS>>();
		
		BufferedReader br = new BufferedReader(new FileReader(file));
		
		while(br.ready()){
			String line = br.readLine();
			String[] split = line.split("\t");
			String gene = split[0];
			String chr = split[2];
			String strand = split[3];
			String txStart = split[4];
			String txEnd = split[5];

			if(!map.containsKey(chr)){
				map.put(chr, new LinkedList<TSS>());
			}
			if(strand.equals("-")){
				map.get(chr).add(new TSS(gene, chr, strand, Integer.parseInt(txEnd)));
			}
			else{
				map.get(chr).add(new TSS(gene, chr, strand, Integer.parseInt(txStart)));
			}
		}
		br.close();

		return map;
	}
	
	private LinkedList<TSS> sortByTSS(LinkedList<TSS> tssinchr){
		TreeMap<Integer, LinkedList<TSS>> sortingstructure = new TreeMap<Integer, LinkedList<TSS>>();
		
		for(Iterator<TSS> it = tssinchr.iterator(); it.hasNext();){
			TSS next = it.next();
			if(!sortingstructure.containsKey(next.tss)){
				sortingstructure.put(next.tss, new LinkedList<TSS>());
			}
			sortingstructure.get(next.tss).add(next);
		}
		
		LinkedList<TSS> rv = new LinkedList<TSS>();
		while(!sortingstructure.isEmpty()){
			rv.addAll(sortingstructure.pollFirstEntry().getValue());
		}
		return rv;
	}
	
	public TSSGene[] getNearestGene(String chr, int start, int end){
		TreeSet<String> nearesttss = new TreeSet<String>();
		int nearesttssd = Integer.MAX_VALUE;


		TSS[] slist = _tsssorted.get(chr);
		if(slist == null){
			return new TSSGene[]{ new TSSGene("null", -1)};
		}
		int min = 0;
		int max = slist.length;
		int mid = max/2;
		while(true){
			mid = min+((max-min)/2);
			TSS midtss = slist[mid];
			int tss = midtss.tss;
			if(tss < start){
				mid = min;
				min = mid+1;
			}
			else if(tss > end){
				mid = max;
				max = mid-1;
			}
			else{
				//Found a match, now scan in both directions until the start position is out of range
				nearesttss.add(midtss.genename);
				nearesttssd = 0;

				for(int i = mid-1; i > -1; i--){
					TSS g = slist[i];
					if(g.tss <= end && g.tss >= start){
						nearesttss.add(g.genename);
					}
					else{
						break;
					}
				}
				for(int i = mid+1; i < slist.length; i++){
					TSS g = slist[i];
					if(g.tss <= end && g.tss >= start){
						nearesttss.add(g.genename);
					}
					else{
						break;
					}
				}
				break;
			}
			
			if(min >= max){
				break;
			}
		}
		
		if(nearesttss.isEmpty()){
			//if no match, look at the nearest tss and add the closest one
			TSS midtss = slist[mid];
			int tss = midtss.tss;
			if(tss < start){ //look right
				if(mid == slist.length-1){
					nearesttss.add(midtss.genename);
					nearesttssd = start-tss;
				}
				else{
					TSS rtss = slist[mid+1];
					int ld = start-tss;
					int rd = rtss.tss-end;
					if(ld < rd){
						nearesttss.add(midtss.genename);
						addSameTSSLeft(nearesttss, slist, mid);
					}
					else if(ld > rd){
						nearesttss.add(rtss.genename);
						addSameTSSRight(nearesttss, slist, mid+1);
					}
					else{
						nearesttss.add(midtss.genename);
						addSameTSSLeft(nearesttss, slist, mid);
						nearesttss.add(rtss.genename);
						addSameTSSRight(nearesttss, slist, mid+1);
					}
					nearesttssd = Math.min(ld, rd);
				}
			}
			else { //look left
				if(mid == 0){
					nearesttss.add(midtss.genename);
					nearesttssd = tss-end;
				}
				else{
					TSS ltss = slist[mid-1];
					int rd = tss-end;
					int ld = start-ltss.tss;
					if(ld < rd){
						nearesttss.add(ltss.genename);
						addSameTSSLeft(nearesttss, slist, mid-1);
					}
					else if(ld > rd){
						nearesttss.add(midtss.genename);
						addSameTSSRight(nearesttss, slist, mid);
					}
					else{
						nearesttss.add(ltss.genename);
						addSameTSSLeft(nearesttss, slist, mid-1);
						nearesttss.add(midtss.genename);
						addSameTSSRight(nearesttss, slist, mid);
					}
					nearesttssd = Math.min(ld, rd);
				}
			}
		}
		
		String[] genes = nearesttss.toArray(new String[0]);
		
		TSSGene[] rv = new TSSGene[genes.length];
		for(int i = 0; i < genes.length; i++){
			rv[i] = new TSSGene(genes[i], nearesttssd);
		}
		
		return rv;
	}
	
	private void addSameTSSLeft(TreeSet<String> set, TSS[] alltss, int position){
		int itss = alltss[position].tss;
		for(int i = position; i > -1; i--){
			int ntss = alltss[i].tss;
			if(itss == ntss){
				set.add(alltss[i].genename);
			}
			else {
				break;
			}
		}
	}
	
	private void addSameTSSRight(TreeSet<String> set, TSS[] alltss, int position){
		int itss = alltss[position].tss;
		for(int i = position; i < alltss.length; i++){
			int ntss = alltss[i].tss;
			if(itss == ntss){
				set.add(alltss[i].genename);
			}
			else{
				break;
			}
		}
	}
	
	public class TSS {
		private String genename;
		private String chr;
		private int tss;
		private String strand;
		
		public TSS(String g, String c, String str, int t){
			genename = g;
			chr = c;
			tss = t;
			strand = str;
		}
		
		public String getGeneName(){
			return genename;
		}
		
		public String getChr(){
			return chr;
		}
		
		public int getTSS(){
			return tss;
		}
		
		public String getStrand(){
			return strand;
		}
	}
}
