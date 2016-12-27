package prioritizationannotation;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map.Entry;
import java.util.TreeMap;
import java.util.TreeSet;

import nearestgene.NearestTSSUtil;
import nearestgene.NearestTSSUtil.TSS;
import nearestgene.TSSGene;


public class Prioritization {
	
	//P - Promoter
	//I - Insulator
	//R - Repressed
	//A - Asleep (Poised) Promoter
	//T - Transcribed/Elongation
	//E - Enhancer
	//S - Signal Low/Repetitive Elements
	
	//GM
	//	0,176,80 -> transcribed  - T
	//	10,190,254 -> insulator - I
	//	127,127,127 -> repressed - R
	//	207,11,198 -> poised promoter - A 
	//	250,202,0 -> strong enhancer - E
	//	255,0,0 -> active promoter - P
	//	255,105,105 -> weak promoter - P
	//	255,252,4 -> weak enhancer - E
	//	255,255,255 -> low signal  - S

	//k562
	//	0,176,80 -> transcribed - T
	//	10,190,254 -> insulator - I
	//	127,127,127 -> repressed - R
	//	207,11,198 -> poised promoter - A
	//	250,202,0 -> strong enhancer - E
	//	255,0,0 -> active promoter - P
	//	255,105,105 -> weak promoter - P
	//	255,252,4 -> weak enhancer - E
	//	255,255,255 -> low signal - S
	
	
	//MCF7
		
	//	1	Elongation - T
	//	2	Enhancer - E
	//	3	Weak_enhancer - E
	//	4	Strong_enhancer - E
	//	5	Strong_enhancer - E
	//	6	Strong_promoter - P
	//	7	Low_signal - S
	//	8	Strong_promoter - P
	//	9	Strong_enhancer - E
	//	10	Repetitive_elements - S
	//	11	Insulator - I
	//	12	Polycomb_repressed - R
	private static int uptss = 2000;
	private static int downtss = 2000;
	
	public static void main(String[] args){
		if(args.length == 9){
			Prioritization p = new Prioritization();
			//args[0] - Filepath to DNase/Node Annotations
			//args[1] - Filepath to ChromHMM States preprocessed where first column is chromosome (eg. chr1), second and third are start and end positions, and the fourth column has a character among { P, I, R, A, T, E, S }
			//args[2] - Filepath to Broad Domains (BED Format/chr star end)
			//args[3] - Filepath to Stretch Enhancers (BED Format/chr star end)
			//args[4] - Filepath to Super Enhancers (BED Format/chr star end)
			p.getAnnotations(args[0], args[1], args[2], args[3], args[4], args[8], args[5]);
			uptss = Integer.parseInt(args[6]);	//TODO There may be issues if the upstream and downstream are different, may need to look into this further
			downtss = Integer.parseInt(args[7]);
			
		}
		else{
			System.out.println("This program requires exactly 9 arguments:");
			System.out.println("1st Argument - Filepath to Sites to Annotate (BED Format/chr star end)");
			System.out.println("2nd Arugment - Filepath to ChromHMM States preprocessed where first column is chromosome (eg. chr1), second and third are start and end positions, and the fourth column has a character among { P, I, R, A, T, E, S } (see below)");
			System.out.println("3rd Argument - Filepath to Broad Domains (BED Format/chr star end)");
			System.out.println("4th Argument - Filepath to Stretch Enhancers (BED Format/chr star end)");
			System.out.println("5th Argument - Filepath to Super Enhancers (BED Format/chr star end)");
			System.out.println("6th Argument - Filepath to Unzipped refflat gene annotations file from UCSC annotation database ftp downloads (http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz)");
			System.out.println("7th Argument - Distance upstream (strand dependent) of a TSS specified as a promoter");
			System.out.println("8th Argument - Distance downstream (strand dependent) of a TSS specified as a promoter");
			System.out.println("9th Argument - Output Directory");

			
			System.out.println("");
			System.out.println("P - Promoter");
			System.out.println("I - Insulator");
			System.out.println("R - Repressed");
			System.out.println("A - Asleep (Poised) Promoter");
			System.out.println("T - Transcribed/Elongation");
			System.out.println("E - Enhancer");
			System.out.println("S - Signal Low/Repetitive Elements");
		}
	}
	
	private BDStats _bdstats;
	
	public void getAnnotations(String dnasefile, String chromhmmfile, String bdfile, String stefile, String sefile, String outdir, String refflatfile){
		
		ChromHMMLocation[] ca = getChromHMMAnnotated(refflatfile, dnasefile, chromhmmfile);
		Location[] bds, stes, ses;
		try {
			bds = readLocations(bdfile);
			stes = readLocations(stefile);
			
			TreeMap<Integer, Integer> allbd = new TreeMap<Integer, Integer>();
			TreeMap<Integer, Integer> abd = new TreeMap<Integer, Integer>();
			_bdstats = new BDStats();
			annotateBD(refflatfile, bds, ca, allbd, abd);
			writeStats(outdir+"_bd_all_stats.txt", allbd);
			writeStats(outdir+"_bd_app_stats.txt", abd);
			writeAnnotationFile(outdir+"_BD.txt", ca, "BD");
			_bdstats.writeFile(outdir+"_single_bd_stats.txt");

			TreeMap<Integer, Integer> allste = new TreeMap<Integer, Integer>();
			TreeMap<Integer, Integer> aste = new TreeMap<Integer, Integer>();
			annotateStE(stes, ca, allste, aste);
			writeStats(outdir+"_ste_all_stats.txt", allste);
			writeStats(outdir+"_ste_app_stats.txt", aste);
			writeAnnotationFile(outdir+"_StE.txt", ca, "StE");

			
			ses = readLocations(sefile);
			TreeMap<Integer, Integer> allse = new TreeMap<Integer, Integer>();
			TreeMap<Integer, Integer> ase = new TreeMap<Integer, Integer>();
			annotateSE(ses, ca, allse, ase);
			writeStats(outdir+"_se_all_stats.txt", allse);
			writeStats(outdir+"_se_app_stats.txt", ase);
			writeAnnotationFile(outdir+"_SE.txt", ca, "SE");

			writeAnnotationFile(outdir+"_P.txt", ca, "P");
			writeAnnotationFile(outdir+"_E.txt", ca, "E");
			writeAnnotationFile(outdir+"_I.txt", ca, "I");
			writeAnnotationFile(outdir+"_PP.txt", ca, "PP");
			writeAnnotationFile(outdir+"_R.txt", ca, "R");
			writeAnnotationFile(outdir+"_T.txt", ca, "T");
			writeAnnotationFile(outdir+"_LS.txt", ca, "LS");

		} catch (Exception e1) {
			e1.printStackTrace();
		}
		
	}
	
	private void writeAnnotationFile(String file, ChromHMMLocation[] dnase, String state) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		for(int i = 0; i < dnase.length; i++){
			ChromHMMLocation cur = dnase[i];
			if(cur.getState().equals(state)){
				bw.write(cur.getChr()+"\t"+cur.getStart()+"\t"+cur.getEnd()+"\n");
			}
		}
		bw.flush();
		bw.close();
	}
	
	private void writeStats(String file, TreeMap<Integer,Integer> stats) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(file));
		while(!stats.isEmpty()){
			Entry<Integer, Integer> e = stats.pollFirstEntry();
			bw.write(e.getKey()+"\t"+e.getValue()+"\n");
		}
		bw.flush();
		bw.close();
	}
	
	//DNASE File - 3 column tab delimited (Chr, Start, End)
	//ChromHMM File - 4 Column Tabe delimited (Chr, Start, End, State) 4th column should be pre-processed with PIRATES
	private ChromHMMLocation[] getChromHMMAnnotated(String refflatfile, String dnasefile, String chromhmmfile){
		
		try {
			Location[] dnase = readLocations(dnasefile);
			ChromHMMLocation[] chromhmm = readChromHMMLocationLocations(chromhmmfile);
			Location[] promoters = getTSS(refflatfile);
			
			return getAnnotatedDNASEPeaks(dnase, chromhmm, promoters);
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}
	
	private void annotateSE(Location[] l, ChromHMMLocation[] dnase, TreeMap<Integer, Integer> allstats, TreeMap<Integer, Integer> astats) throws Exception{
		TreeMap<String, Location[]> dnasesorted = getChrSortedLocations(dnase);
		for(int i = 0; i < l.length; i++){
			annotateSE(l[i], dnasesorted.get(l[i].getChr()), allstats, astats);
		}
	}
	
	private void annotateSE(Location l, Location[] dnase, TreeMap<Integer, Integer> allstats, TreeMap<Integer, Integer> astats) throws Exception{
		int start = l.getStart();
		int end = l.getEnd();
		
		int all = 0;
		LinkedList<ChromHMMLocation> candidates = new LinkedList<ChromHMMLocation>();
		if(dnase != null){
			int fi = findFirstIndex(dnase, l);
			if(fi != -1){
				for(int i = fi; i < dnase.length; i++){
					ChromHMMLocation cur = (ChromHMMLocation) dnase[i];
					int cs = cur.getStart();
					int ce = cur.getEnd();
					if(ce < start){
						continue;
					}
					if(end <= cs){
						break;
					}
					
					all++;
					if(cur.getState().equals("E") || cur.getState().equals("StE")){
						candidates.add(cur);
					}
				}
			}
		}
		
		if(!allstats.containsKey(all)){
			allstats.put(all, 0);
		}		
		allstats.put(all, allstats.get(all)+1);

		int numca = candidates.size();
		if(!astats.containsKey(numca)){
			astats.put(numca, 0);
		}
		astats.put(numca, astats.get(numca)+1);
		
		ChromHMMLocation[] ca = candidates.toArray(new ChromHMMLocation[0]);
		for(int i = 0; i < ca.length; i++){
			ca[i].setState("SE");
		}

	}
	
	private void annotateStE(Location[] l, ChromHMMLocation[] dnase, TreeMap<Integer, Integer> allstats, TreeMap<Integer, Integer> astats) throws Exception{
		TreeMap<String, Location[]> dnasesorted = getChrSortedLocations(dnase);
		for(int i = 0; i < l.length; i++){
			annotateStE(l[i], dnasesorted.get(l[i].getChr()), allstats, astats);
		}
	}
	
	private void annotateStE(Location l, Location[] dnase,  TreeMap<Integer, Integer> allstats, TreeMap<Integer, Integer> astats) throws Exception{
		int start = l.getStart();
		int end = l.getEnd();
		
		int all = 0;
		LinkedList<ChromHMMLocation> candidates = new LinkedList<ChromHMMLocation>();
		if(dnase != null){
			int fi = findFirstIndex(dnase, l);
			if(fi != -1){
				for(int i = fi; i < dnase.length; i++){
					ChromHMMLocation cur = (ChromHMMLocation) dnase[i];
					int cs = cur.getStart();
					int ce = cur.getEnd();
					if(ce < start){
						continue;
					}
					if(end <= cs){
						break;
					}
					
					all++;
					if(cur.getState().equals("E")){
						candidates.add(cur);
					}
				}
			}
		}
		
		if(!allstats.containsKey(all)){
			allstats.put(all, 0);
		}		
		allstats.put(all, allstats.get(all)+1);

		int numca = candidates.size();
		if(!astats.containsKey(numca)){
			astats.put(numca, 0);
		}
		astats.put(numca, astats.get(numca)+1);
		
		ChromHMMLocation[] ca = candidates.toArray(new ChromHMMLocation[0]);
		for(int i = 0; i < ca.length; i++){
			ca[i].setState("StE");
		}
		

	}
	
	
	private void annotateBD(String refflatfile, Location[] l, ChromHMMLocation[] dnase, TreeMap<Integer, Integer> allstats, TreeMap<Integer, Integer> astats) throws Exception{
		TreeMap<String, Location[]> dnasesorted = getChrSortedLocations(dnase);
		TreeMap<String, Location[]> promotersorted = getChrSortedLocations(getTSS(refflatfile));

		NearestTSSUtil ntss = new NearestTSSUtil(refflatfile);

		for(int i = 0; i < l.length; i++){
			String chr = l[i].getChr();
			annotateBD(l[i], dnasesorted.get(chr), allstats, astats, promotersorted.get(chr), ntss);
		}
	}
	
	private void annotateBD(Location l, Location[] dnase,  TreeMap<Integer, Integer> allstats, TreeMap<Integer, Integer> astats, Location[] promoters, NearestTSSUtil ntss) throws Exception{
		int start = l.getStart();
		int end = l.getEnd();
		int d = end-start;
		
		int all = 0;
		LinkedList<ChromHMMLocation> candidates = new LinkedList<ChromHMMLocation>();
		if(dnase != null){
			int fi = findFirstIndex(dnase, l);
			if(fi != -1){
				for(int i = fi; i < dnase.length; i++){
					ChromHMMLocation cur = (ChromHMMLocation) dnase[i];
					int cs = cur.getStart();
					int ce = cur.getEnd();
					if(ce < start){
						continue;
					}
					if(end <= cs){
						break;
					}
					
					all++;
					if(cur.getState().equals("P")){
						candidates.add(cur);
					}
				}
			}
		}
		
		if(!allstats.containsKey(all)){
			allstats.put(all, 0);
		}		
		allstats.put(all, allstats.get(all)+1);

		int numca = candidates.size();
		if(!astats.containsKey(numca)){
			astats.put(numca, 0);
		}
		astats.put(numca, astats.get(numca)+1);
		
		
		ChromHMMLocation[] ca = candidates.toArray(new ChromHMMLocation[0]);
		int[] overlaps = new int[ca.length];
		int[] dists = new int[ca.length];
		int maxoverlap = 0;
		int mindist = Integer.MAX_VALUE;
		for(int i = 0; i < ca.length; i++){
			ChromHMMLocation cur = ca[i];
			int cstart = cur.getStart();
			int cend = cur.getEnd();
			int overlap;
			if(cstart > start && cend < end){
				overlap = cend-cstart;
			}
			else if(cstart > start){
				overlap = end-cstart;
			}
			else{
				overlap = cend-start;
			}
			
			TSSGene[] genes = ntss.getNearestGene(cur.getChr(), cstart, cend);
			dists[i] = genes[0].getDistance();
			
			maxoverlap = Math.max(maxoverlap, overlap);
			mindist = Math.min(mindist, dists[i]);
			overlaps[i] = overlap;
		}
		
		LinkedList<ChromHMMLocation> maxoverlapl = new LinkedList<ChromHMMLocation>();
		for(int i = 0; i < overlaps.length; i++){
			if(overlaps[i] == maxoverlap){
				if(dists[i] == mindist){
					maxoverlapl.add(ca[i]);
				}
			}
		}
		
		if(ca.length > 0){
			if(maxoverlapl.size() == 0){
				_bdstats.miss++;
				int smaxoverlap = 0;
				for(int i = 0; i < overlaps.length; i++){
						if(dists[i] == mindist){
							smaxoverlap = Math.max(smaxoverlap, overlaps[i]);
						}
				}
				for(int i = 0; i < overlaps.length; i++){
					if(dists[i] == mindist && smaxoverlap == overlaps[i]){
						maxoverlapl.add(ca[i]);
						_bdstats.choicestats.add(new int[]{dists[i], overlaps[i], d});
					}
				}
	
			}
			else{
				if(ca.length > 1){
					_bdstats.hit++;
				}
				else{
					_bdstats.nochoice++;
				}
			}
		}
		else{
			_bdstats.notin++;
		}
		
		if(maxoverlapl.size() > 1){
			//Need to see which one overlaps the most? 
			//throw new Exception("That happened."); //Only going to code this if it actually occurs
			System.out.println("BD:"+maxoverlapl.size()+" undecided, annotating all|"+l.getChr()+":"+start+"-"+end);
			for(Iterator<ChromHMMLocation> it = maxoverlapl.iterator(); it.hasNext();){
				it.next().setState("BD");
			}
		}
		else{
			if(!maxoverlapl.isEmpty()){
				maxoverlapl.getFirst().setState("BD");
			}
		}
	}
	
	private ChromHMMLocation[] getAnnotatedDNASEPeaks(Location[] dnase, ChromHMMLocation[] chromhmm, Location[] promoters) throws Exception{
		TreeMap<String, Location[]> chromhmmsorted = getChrSortedLocations(chromhmm);
		TreeMap<String, Location[]> promotersorted = getChrSortedLocations(promoters);
		
		LinkedList<ChromHMMLocation> rv = new LinkedList<ChromHMMLocation>();
		for(int i = 0; i < dnase.length; i++){
			String chr = dnase[i].getChr();
			if(chromhmmsorted.containsKey(chr)){
				rv.add(getAnnotatedDNASEPeak(dnase[i], chromhmmsorted.get(chr), promotersorted.get(chr)));
			}
		}
		return rv.toArray(new ChromHMMLocation[0]);
	}
	
	private ChromHMMLocation getAnnotatedDNASEPeak(Location l, Location[] chromhmm, Location[] promoters){
		String chr = l.getChr();
		int start = l.getStart();
		int end = l.getEnd();
		
		boolean p = false, e = false, ins = false, pp = false, r = false, t = false, ls = false;
		
		int fi = findFirstIndex(chromhmm, l);
		if(fi == -1){
			System.out.println("Error: ChromHMM Location Missing: "+chr+":"+start+"-"+end);
		}
		else{
			for(int i = fi; i < chromhmm.length; i++){
				ChromHMMLocation cur = (ChromHMMLocation) chromhmm[i];
				int cs = cur.getStart();
				int ce = cur.getEnd();
				if(ce < start){
					continue;
				}
				if(end <= cs){
					break;
				}
				
				if(cur.getState().equals("P")){
					p = true;
				}
				else if(cur.getState().equals("I")){
					ins = true;
				}
				else if(cur.getState().equals("R")){
					r = true;
				}
				else if(cur.getState().equals("A")){
					pp = true;
				}
				else if(cur.getState().equals("T")){
					t = true;
				}
				else if(cur.getState().equals("E")){
					e = true;
				}
				else if(cur.getState().equals("S")){
					ls = true;
				}
			}
		}
		
		return getAnnotation(l, p, e, ins, pp, r, t, ls, promoters);
	}
	
	private ChromHMMLocation getAnnotation(Location l, boolean p, boolean e, boolean i, boolean pp, boolean r, boolean t, boolean ls, Location[] promoters){
		ChromHMMLocation rv = new ChromHMMLocation(-1, l.getChr(), l.getStart(), l.getEnd(), "");
		if(p && e){
			if(findFirstIndex(promoters, l) == -1){
				rv.setState("E");
			}
			else{
				rv.setState("P");
			}
		}
		else{
			if(p){
				rv.setState("P");
			}
			else if(e){
				rv.setState("E");
			}
			else if(i){
				rv.setState("I");
			}
			else if(pp){
				rv.setState("PP");
			}
			else if(r){
				rv.setState("R");
			}
			else if(t){
				rv.setState("T");
			}
			else {
				rv.setState("LS");
			}
		}
		return rv;
	}
	
	private Location[] getTSS(String refflatpath){
		NearestTSSUtil ntssu = null;
		try {
			ntssu = new NearestTSSUtil(refflatpath);
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		TSS[] tsslocations = ntssu.getTSSLocations();
		
		TreeSet<String> seen = new TreeSet<String>();
		LinkedList<Location> rv = new LinkedList<Location>();
		for(int i = 0; i < tsslocations.length; i++){
			String chr = tsslocations[i].getChr();
			int tss = tsslocations[i].getTSS();
			if(tsslocations[i].getStrand().equals("-")){
				String seenkey = chr+(tss-downtss);
				if(!seen.contains(seenkey)){
					rv.add(new Location(-1, chr, tss-downtss, tss+uptss));
					seen.add(seenkey);
				}
			}
			else{
				String seenkey = chr+(tss-uptss);
				if(!seen.contains(seenkey)){
					rv.add(new Location(-1, chr, tss-uptss, tss+downtss));
					seen.add(seenkey);
				}
			}
		}
		
		return rv.toArray(new Location[0]);
	}
	
	
	
	//Returns -1 if there is no position that overlaps the region
	private int findFirstIndex(Location[] loci, Location l){
		int start = l.getStart();
		int end = l.getEnd();
		
		//Binary search to get an approximate index
		int s = 0;
		int e = loci.length-1;
		while((e-s) > 1){
			int mi = s+((e-s)/2);
			Location ml = loci[mi];
			int mstart = ml.getStart();
			if(mstart < start){
				s = mi;
			}
			else if(mstart > start){
				e = mi;
			}
			else{
				s = mi;
				e = mi;
				break;
			}
		}
		
		//Loop to make sure we have the location just before the start position goes beyond
		for(int i = s+1; i < loci.length; i++){
			if(loci[i].getStart() > start){
				s = i-1;
				break;
			}
		}
		
		//Choose position just before or just after depending on the overlap. 
		//If neither of them do return -1
		int n = s+1;
		if(s < loci.length && loci[s].getEnd() >= start){
			return s;
		}
		else if(n < loci.length && loci[n].getStart() <= end && loci[n].getEnd() >= start){
			return n;
		}
		else {
			return -1;
		}
	}
	
	private Location[] readLocations(String file) throws IOException{
		LinkedList<Location> rv = new LinkedList<Location>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		int id = 0;
		while(br.ready()){
			String[] split = br.readLine().split("\t");
			if(split.length > 2){
				rv.add(new Location(id++, split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2])));
			}
		}
		br.close();
		return rv.toArray(new Location[0]);
	}
	
	private ChromHMMLocation[] readChromHMMLocationLocations(String file) throws IOException{
		LinkedList<ChromHMMLocation> rv = new LinkedList<ChromHMMLocation>();
		BufferedReader br = new BufferedReader(new FileReader(file));
		while(br.ready()){
			String[] split = br.readLine().split("\t");
			if(split.length > 4){
				rv.add(new ChromHMMLocation(-1, split[0], Integer.parseInt(split[1]), Integer.parseInt(split[2]), split[3]));
			}
		}
		br.close();
		return rv.toArray(new ChromHMMLocation[0]);
	}
	
	private TreeMap<String, Location[]> getChrSortedLocations(Location[] loci) throws Exception{
		TreeMap<String, LinkedList<Location>> m = new TreeMap<String, LinkedList<Location>>();
		for(int i = 0; i < loci.length; i++){
			String key = loci[i].getChr();
			if(!m.containsKey(key)){
				m.put(key, new LinkedList<Location>());
			}
			m.get(key).add(loci[i]);
		}
		TreeMap<String, Location[]> rv = new TreeMap<String, Location[]>();
		while(!m.isEmpty()){
			Entry<String, LinkedList<Location>> e = m.pollFirstEntry();
			rv.put(e.getKey(), getStartSortedLocations(e.getValue().toArray(new Location[0])));
		}
		return rv;
	}
	
	private Location[] getStartSortedLocations(Location[] loci) throws Exception{
		TreeMap<Integer, Location> m = new TreeMap<Integer, Location>();
		for(int i = 0; i < loci.length; i++){
			int key = loci[i].getStart();
			if(m.containsKey(key)){
				throw new Exception("Duplicate Start Position!"); //Assumes that there is not a location that has the same start position
			}
			else{
				m.put(key, loci[i]);
			}
		}
		return m.values().toArray(new Location[0]);
	}

	private class BDStats {
		public int notin = 0;
		public int nochoice = 0;
		public int hit = 0;
		public int miss = 0;
		public LinkedList<int[]> choicestats;
		public BDStats() {
			choicestats = new LinkedList<int[]>();
		}
		
		public void writeFile(String file) throws IOException{
			BufferedWriter bw = new BufferedWriter(new FileWriter(file));
			bw.write("Not in:\t"+notin+"\n");
			bw.write("No Choice:\t"+nochoice+"\n");
			bw.write("Max Overlap Matches Min TSS Dist:\t"+hit+"\n");
			bw.write("Misses:\t"+miss+"\n");
			for(Iterator<int[]> it = choicestats.iterator(); it.hasNext();){
				int[] next = it.next();
				bw.write(next[0]+"\t"+next[1]+"\t"+next[2]+"\n");
			}
			bw.flush();
			bw.close();
		}
	}
	
}
