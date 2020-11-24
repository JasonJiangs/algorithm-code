package as8;

import java.util.Arrays;

public class Knapsack { // knapsack problem

	public static void main(String[] args) {
		long totalTime = 0;
		for (int n = 0; n <= 150; n = n + 1) { // gap is 1, every time increase 1, from 0 to 150
			int[] weights = new int[n];
			int[] profits = new int[n];
			totalTime = 0;
			for (int count1 = 0; count1 < 50; count1++) { // repeat 50 times

				// Create the random number
				for (int index = 0; index < n; index++) {
					weights[index] = (int) (Math.random() * 50);
					profits[index] = (int) (Math.random() * 100);
				}
				// Analysis the total time
				long startTime = System.nanoTime();
				dynamic_Programming(weights, profits, 60);
				//bruteForce_Recursion(weights.length-1, 60, weights, profits);
				totalTime = totalTime + System.nanoTime() - startTime;
			}
			System.out.println(totalTime / 50);
		}
	}

	// dynamic programming to solve the problem
	public static int dynamic_Programming(int[] weight, int[] value, int capacity) {
		int weightLen = weight.length;
		int valueLen = capacity + 1;
		int maxValue = 0;
		int[][] v = new int[weightLen][valueLen];
		for (int i = 0; i < weightLen; i++) {
			for (int j = 0; j < valueLen; j++) {
				if (i == 0) {
					v[i][j] = 0;
				} else if (j == 0) {
					v[i][j] = 0;
				} else {
					if (weight[i] > j) {
						v[i][j] = v[i - 1][j];
					} else if (weight[i] <= j) {
						v[i][j] = Math.max(v[i - 1][j], v[i - 1][j - weight[i]] + value[i]);
					}
					maxValue = v[i][j];
				}
			}
		}
		return maxValue;
	}

	// Recursion form, i means i-th item, j means j space left.
	public static int bruteForce_Recursion(int i, int j, int[] weights, int[] values) {
		int r = 0;
		if (i == -1) {
			return 0;
		}
		// If left space is larger than the items we can put,
		// here we would have two choices, put or not put,
		// then recursively split the question.
		if (j >= weights[i]) {
			int r1 = bruteForce_Recursion(i - 1, j - weights[i], weights, values) + values[i]; // put i-th item
			int r2 = bruteForce_Recursion(i - 1, j, weights, values);// don't put i-th item
			// Compare these two condition's result and get the larger value under this
			// layer.
			r = Math.max(r1, r2);
		}
		return r;
		// r will be returned for recursive purpose or the final result.
	}
}
