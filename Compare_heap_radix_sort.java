package mypackage;

public class Compare_heap_radix_sort {
		public static void main(String[] args) {
//			// create array manually
//			int[] size_arr = { 10, 50, 100, 500, 1000, 5000, 10000, 15000, 20000, 25000, 30000, 35000, 40000, 45000,
//					50000 };
			int[] size_arr = new int[500];
			for (int i = 0; i < size_arr.length; i++) {
				size_arr[i] = (i + 1) * 100;
			}
			// execute comparison function and print out the result
			for (int i = 0; i < size_arr.length; i++) {
				comparison(size_arr[i]);
			}
		}

		// the main method, create random array, repeatedly execute and calculate cost,
		// then average them
		public static void comparison(int Size) {
			int k = 1000000000; // maximum number k
			int repetition = 100; // times of executions
			long[] time_heapSort = new long[repetition]; // calculation of mean for the heap sort
			long[] time_radixSort = new long[repetition]; // calculation of mean for the radix sort

			for (int i = 0; i < repetition; i++) { // repeatedly get time cost array
				// copy array for radix sort to use
				int[] arr = randomArr(Size, k);
				int[] temp = copyArr(arr);
				// calculate the time cost of heap sort, single result stored
				long startTime1 = System.nanoTime();
				heapSort(arr);
				long endTime1 = System.nanoTime();
				long time_taken1 = endTime1 - startTime1;
				time_heapSort[i] = time_taken1;
				// calculate the time cost of radix sort, single result stored
				long startTime2 = System.nanoTime();
				radixSort(temp);
				long endTime2 = System.nanoTime();
				long time_taken2 = endTime2 - startTime2;
				time_radixSort[i] = time_taken2;
			}
			// calculate the average time cost for each sort function
			long averageT_HeapSort = doubleArrAverage(time_heapSort);
			long averageT_RadixSort = doubleArrAverage(time_radixSort);
			//System.out.println(averageT_RadixSort);
			System.out.println(averageT_HeapSort);
		}

		// copy the created random array for the other sort method to sort
		public static int[] copyArr(int arr[]) {
			int copy_array[] = new int[arr.length];
			copy_array = arr;
			return copy_array;
		}

		// when repeat several times and get the results array, calculate the mean
		// values
		public static long doubleArrAverage(long[] arr) {
			long sum = 0;
			for (int i = 0; i < arr.length; i++) {
				sum += arr[i];
			}
			return sum / arr.length;
		}

		// method for creating a random array
		public static int[] randomArr(int size, int max) {
			// size is the size of random, max is the maximum number which is 100.
			int[] arr = new int[size];
			for (int i = 0; i < arr.length; i++) {
				int x = (int) (Math.random() * max);
				arr[i] = x;
			}
			return arr;
		}

		public static void heapSort(int[] arr) {
			// start from non-child node, from bottom to top, from right to left to adjust
			for (int i = (arr.length - 1) / 2; i >= 0; i--) { // build heap
				adjustHeap(arr, i, arr.length);
			}
			// Adjusting heap structure and swapping heap top and bottom elements
			for (int i = arr.length - 1; i > 0; i--) {
				int temp = arr[i]; // Swap the top element of the heap with the last element of the heap.
				arr[i] = arr[0];
				arr[0] = temp;

				adjustHeap(arr, 0, i);// Re-adjustment of heap
			}
		}

		public static void adjustHeap(int[] arr, int parent, int length) {
			int temp = arr[parent];// treat temp as the parent node
			int lChild = 2 * parent + 1;// left child

			while (lChild < length) {
				int rChild = lChild + 1;// right child
				// If there is a right child node and the value of the right child node is
				// greater than the left child node, then the right child node is selected
				if (rChild < length && arr[lChild] < arr[rChild]) {
					lChild++;
				}
				// If the value of the parent node is already greater than the value of the
				// child node,
				// it ends directly
				if (temp >= arr[lChild]) {
					break;
				}
				// Assign the value of the child node to the parent node.
				arr[parent] = arr[lChild];
				// Select the left child node of the child node and continue filtering down.
				parent = lChild;
				lChild = 2 * lChild + 1;
			}
			arr[parent] = temp;
		}

		private static void radixSort(int[] arr) {
			int max = arr[0]; // maximum number waiting for sort
			int exp; // exponent
			for (int anArr : arr) { // calculate the maximum number
				if (anArr > max) {
					max = anArr;
				}
			}
			for (exp = 1; max / exp > 0; exp *= 10) { // sort start from units digits
				int[] temp = new int[arr.length]; // temp array for waiting sorted elements
				int[] buckets = new int[10]; // ten buckets from 0 to 9

				for (int value : arr) { // store the appearance times of the elements into the buckets
					buckets[(value / exp) % 10]++;// (value / exp) % 10 : the last digit of the value
				}
				for (int i = 1; i < 10; i++) {// change buckets[i]
					buckets[i] += buckets[i - 1];
				}
				for (int i = arr.length - 1; i >= 0; i--) {// store elements into the temp array
					temp[buckets[(arr[i] / exp) % 10] - 1] = arr[i];
					buckets[(arr[i] / exp) % 10]--;
				}
				System.arraycopy(temp, 0, arr, 0, arr.length);// assign the sorted elements to arr
			}
		}
	
}
