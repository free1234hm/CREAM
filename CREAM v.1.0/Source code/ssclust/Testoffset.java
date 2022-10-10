package ssclust;

public class Testoffset {
	public static void main(String[] args) {
		new Testoffset();
	}
	public Testoffset() {
		Gettrc b = new Gettrc();
		double[][] weight = {{1, 1,1, 1,1, 1,1, 1},
		                     {1, 1,1, 1,1, 1,1, 1},
		                     {1, 1,1, 1,1, 1,1, 1},
		                     {1, 1,1, 1,1, 1,1, 1},
		                     {1, 1,1, 1,0, 0,0, 0},
		                     {1, 1,1, 1,0, 0,0, 0},
		                     {1, 1,1, 1,0, 0,0, 0}};

		double trc = b.gettrc(weight);
		System.out.println(trc);
	}
}
