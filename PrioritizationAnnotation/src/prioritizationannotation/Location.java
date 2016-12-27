package prioritizationannotation;

public class Location {
	private int _id;
	private String _chr;
	private int _start;
	private int _end;
	
	public Location(int id, String chr, int start, int end){
		_id = id;
		_chr = chr;
		_start = start;
		_end = end;
	}
	
	public int getId(){
		return _id;
	}
	
	public String getChr(){
		return _chr;
	}
	
	public int getStart(){
		return _start;
	}
	
	public int getEnd(){
		return _end;
	}
	
	public void setChr(String chr){
		_chr = chr;
	}
	
	public void setStart(int s){
		_start = s;
	}
	
	public void setEnd(int e){
		_end = e;
	}
}
