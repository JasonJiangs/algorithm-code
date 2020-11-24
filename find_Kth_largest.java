package mypackage;

import java.util.Random;
import mypackage.test;

public class find_Kth_largest {
	public static void main(String[] args) {
		int[] size_arr = new int[1000]; //500
		for (int i = 0; i < size_arr.length; i++) {
			size_arr[i] = (i + 1) * 100;
		}
		// execute comparison function and print out the result
		for (int i = 0; i < size_arr.length; i++) {
			int i_th = (int) (Math.random() * size_arr[i]);
			comparison(size_arr[i], i_th);
		}
	}

///////////////////////////////////////////////////////////////////////////////////////
	public static void comparison(int Size, int i_th) {
		int k = 1000; // maximum number k
		int repetition = 50;
		long[] time_RandomSelect = new long[repetition]; // calculation of mean for the heap sort
		long[] time_MoMselect = new long[repetition]; // calculation of mean for the radix sort
		for (int i = 0; i < repetition; i++) { // repeatedly get time cost array
			// copy array for selection with median of median to use
			int[] arr = randomArr(Size, k);
			int[] temp = copyArr(arr);
			// calculate the time cost of random select, single result stored
			long startTime1 = System.nanoTime();
			Quick_Select.findKthLargest(arr, i_th);
			long endTime1 = System.nanoTime();
			long time_taken1 = endTime1 - startTime1;
			time_RandomSelect[i] = time_taken1;
			// calculate the time cost of selection with median of median, single result
			// stored
			long startTime2 = System.nanoTime();
			MOM.MedianOfMedian(temp, 0, temp.length-1, i_th);
			long endTime2 = System.nanoTime();
			long time_taken2 = endTime2 - startTime2;
			time_MoMselect[i] = time_taken2;
		}
		// calculate the average time cost for each sort function
		long averageT_RandomSelect = doubleArrAverage(time_RandomSelect);
		long averageT_MoMselect = doubleArrAverage(time_MoMselect);
		System.out.println(averageT_MoMselect);
		//System.out.println(averageT_RandomSelect);
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
}

class Quick_Select {
    public static void swap(int[]nums, int a, int b) {
        int tmp = nums[a];
        nums[a] = nums[b];
        nums[b] = tmp;
    }
    public static int partition(int[] nums, int left, int right, int pivot_index) {
        int pivot = nums[pivot_index];
        // 1. move pivot to end
        swap(nums, pivot_index, right);//不让pivot比较它自己，先移走，回头再移回来
        int store_index = left;
        // 2. move all smaller elements to the left
        for (int i = left; i <= right; i++) {//遍历left - right
            if (nums[i] < pivot) {
                swap(nums, store_index, i);
                store_index++;//每出现一个比pivot小的，就交换，store_index+1，最后的结果nums[store_index]左边都比pivot小
            }
        }
        // 3. move pivot to its final place
        swap(nums, store_index, right);

        return store_index;
    }
    /**
     * 快速选择算法
     * @param left 数组左侧
     * @param right 数组右侧
     * @param k_smallest k
     * @return 第k小的元素
     */
    public static int quickselect(int[] nums, int left, int right, int k_smallest) {
        /*
        Returns the k-th smallest element of list within left..right.
        */
        if (left == right) // If the list contains only one element,
            return nums[left];  // return that element
        // select a random pivot_index
        Random random_num = new Random();
        int pivot_index = left + random_num.nextInt(right - left);

        pivot_index = partition(nums, left, right, pivot_index);

        // the pivot is on (N - k)th smallest position
        if (k_smallest == pivot_index)
            return nums[k_smallest];
            // go left side
        else if (k_smallest < pivot_index)
            return quickselect(nums, left, pivot_index - 1, k_smallest);
        // go right side
        return quickselect(nums, pivot_index + 1, right, k_smallest);
    }
    /**
     * 选择第k大的元素
     * @param nums
     * @param k
     * @return 第k大的元素
     */
    public static int findKthLargest(int[] nums, int k) {
        int size = nums.length;
        // kth largest is (N - k)th smallest
        return quickselect(nums, 0, size - 1, size - k);
    }

    public static void main(String[] args) {
        QuickSelect quickSelect = new QuickSelect();
        int[] nums = {4,3,2,6,5,1,9,9,7,8};
        int kthLargest = quickSelect.findKthLargest(nums, 11);
        System.out.println(kthLargest);
    }
}

class MOM {
	public static void quickSort(int[] a, int left, int right) {
		int i, j, t, pivot;
		if (left > right) {
			return;
		}
		pivot = a[left];
		i = left; // too big index
		j = right;// too small index
		while (i < j) {
			while (a[i] <= pivot && i < j) {
				i++;
			}
			while (a[j] > pivot) {
				j--;
			}
			if (i < j) {
				t = a[i];
				a[i] = a[j];
				a[j] = t;
			}
		}
		a[left] = a[j];// 结束后要把pivot提到中间来
		a[j] = pivot;
		quickSort(a, left, j - 1);
		quickSort(a, j + 1, right);
	}
	public static int MedianOfMedian(int[] a, int p, int q, int i) {
		int x, r, k;
		if (p == q) {
			return a[p];
		}
		x = GoodPivot(a, p, q);
		r = PartitionWithPivot(a, p, q, x);
		k = r - p + 1;
		if (i == k) {
			return a[r];
		}
		if (i < k) {
			return MedianOfMedian(a, p, r - 1, i);
		} else {
			return MedianOfMedian(a, r + 1, q, i - k);
		}
	}

	public static int PartitionWithPivot(int[] a, int p, int q, int x) {
		int i = p;
		int j = q;
		int t;
		int pivot = x;
		int findindex = 0;
		boolean find = false;
		while (i < j) {
			while (a[i] <= pivot && i < j) {
				if (find == false) {
					if (a[i] == pivot) {
						findindex = i;
						find = true;
					}
				}
				i++;
			}
			while (a[j] > pivot) {
				j--;
			}
			if (i < j) {
				t = a[i];
				a[i] = a[j];
				a[j] = t;
			}
		}
		a[findindex] = a[j];// 结束后要把pivot提到中间来
		a[j] = pivot;
		return j;
	}

	public static int GoodPivot(int[] a, int p, int q) {
		int[] b = new int[(q - p + 5) / 5];
		for (int j = p; j <= q; j += 5) {
			b[(j - p) / 5] = MedianOfArray(a, j, Math.min(j + 4, q));
		}
		return MedianOfArray(b, 0, b.length - 1);
	}

	public static int MedianOfArray(int[] a, int p, int q) {
		quickSort(a, p, q);
		return a[(p + q) / 2];
	}
}