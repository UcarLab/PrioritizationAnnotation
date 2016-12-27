package prioritizationannotation;

public class ChromHMMLocation extends Location {

	private String _state;

	public ChromHMMLocation(int id, String chr, int start, int end, String state) {
		super(0, chr, start, end);
		_state = state;
	}

	public String getState(){
		return _state;
	}
	
	public void setState(String state){
		_state = state;
	}
}
