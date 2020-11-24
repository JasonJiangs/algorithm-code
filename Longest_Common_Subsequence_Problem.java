package as9;

public class  Longest_Common_Subsequence_Problem {

	public static void main(String[] args) {
		// worst running time is O(n*2^m)
		String chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
		int m = 0; // n for arr1, m for arr2
		int n = 2; // fix n, increase m, which can create exponential growth
		String str1 = "";
		for (int i = 0; i < n; i++) {
			str1 = str1 + chars.charAt((int) (Math.random() * 26));
		}
		String str2 = "";
		String str3 = "";
		for (int i = 0; i < 1000; i++) {
			m = m + 1;
			for (int j = 0; j < m; j++) {
				str2 = str2 + chars.charAt((int) (Math.random() * 26));
				str3 = str3 + chars.charAt((int) (Math.random() * 26));
			}
			long startTime = System.nanoTime();
			
//			int result = brutal_force2.brutalForce(str1, str2); // implement brute force solution
			
			int[][] re = buttom_up.longestCommonSubsequence(str2, str3);
			buttom_up.print(re, str2, str3, str2.length(), str3.length());
			int result = buttom_up.strLength;
			
//			int result = top_down.topDown(str2, str3, str2.length()-1, str3.length()-1, 0);
			
//			TopDown d = new TopDown(str2.length() + 1, str3.length() + 1);
	//		int result = d.TopDown(str2, str3, str2.length() - 1, str3.length() - 1);
			long closeTime = System.nanoTime();
			long timeConsume = closeTime - startTime;
		//	System.out.print("For strings " + str3 + " and " + str2 + ", it costs " + timeConsume);
	//		System.out.println(" nanotimes. Result is " + result + ".");
			System.out.println(timeConsume);
			str2 = "";
			str3 = "";
		}
	}

}

class TopDown {
	static int[][] storage;

	public TopDown(int m, int n) {
		this.storage = new int[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				storage[i][j] = -1;
			}
		}
	}

	public static int TopDown(String s1, String s2, int i, int j) {
		if (i < 0 || j < 0) {
			return 0;
		}
		if (s1.charAt(i) == s2.charAt(j)) {
			if (storage[i][j] == -1) {
				storage[i][j] = TopDown(s1, s2, i - 1, j - 1);
			}
			return 1 + storage[i][j];
		} else {
			if (storage[i][j + 1] == -1) {
				storage[i][j + 1] = TopDown(s1, s2, i - 1, j);
			}
			if (storage[i + 1][j] == -1) {
				storage[i + 1][j] = TopDown(s1, s2, i, j - 1);
			}
			return Math.max(storage[i][j + 1], storage[i + 1][j]);
		}
	}
}
