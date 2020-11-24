package as11;

import java.util.Arrays;
import java.util.Random;

public class MST_Algorithm {
	private int mEdgNum; // number of edges
	private String[] mVexs; // all the names of the vertexes
	private int[][] mMatrix; // adjacent matrix
	private static final int INF = Integer.MAX_VALUE; // use maximum value of integer

	// create graph using existed matrix and vertexes
	public MST_Algorithm(String[] vexs, int[][] matrix) {
		// initialize the vertexes
		int vlen = vexs.length;
		mVexs = new String[vlen];
		for (int i = 0; i < mVexs.length; i++) {
			mVexs[i] = vexs[i];
		}
		// initialize edges
		mMatrix = new int[vlen][vlen];
		for (int i = 0; i < vlen; i++) {
			for (int j = 0; j < vlen; j++) {
				mMatrix[i][j] = matrix[i][j];
			}
		}
		// count number of edges
		mEdgNum = 0;
		for (int i = 0; i < vlen; i++) {
			for (int j = i + 1; j < vlen; j++) {
				if (mMatrix[i][j] != INF) {
					mEdgNum++;
				}
			}
		}
	}

	// return the position of the target vertex
	private int getPosition(String aim) {
		for (int i = 0; i < mVexs.length; i++) {
			if (mVexs[i].equals(aim)) {
				return i;
			}
		}
		return -1;
	}

	// return the index of the first adjacent vertex of v, if fail then return -1
	private int firstVertex(int v) {
		if (v < 0 || v > (mVexs.length - 1)) {
			return -1;
		}
		for (int i = 0; i < mVexs.length; i++) {
			if (mMatrix[v][i] != 0 && mMatrix[v][i] != INF) {
				return i;
			}
		}
		return -1;
	}

	// return index of v relative to the next adjacent vertex of w, return -1 if
	// fail
	private int nextVertex(int v, int w) {
		if (v < 0 || v > (mVexs.length - 1) || w < 0 || w > (mVexs.length - 1)) {
			return -1;
		}
		for (int i = w + 1; i < mVexs.length; i++) {
			if (mMatrix[v][i] != 0 && mMatrix[v][i] != INF) {
				return i;
			}
		}
		return -1;
	}

	// print the matrix graph in a strict format
	public void print() {
		System.out.printf("Martix Graph:\n");
		for (int i = 0; i < mVexs.length; i++) {
			for (int j = 0; j < mVexs.length; j++) {
				System.out.printf("%10d ", mMatrix[i][j]);
			}
			System.out.printf("\n");
		}
	}

	// MST using prim algorithm (array)
	// start from the index of vertex of center
	public void prim_Array(int center) {
		int num = mVexs.length; // number of vertexes
		int index = 0; // index of the MST
		String[] prims = new String[num]; // result array for MST
		int[] weights = new int[num]; // weights between vertexes
		// assign the first vertex with the center vertex
		prims[index++] = mVexs[center];
		// initialize the weight array for vertexes
		// store weights from every vertex to the center vertex
		for (int i = 0; i < num; i++) {
			weights[i] = mMatrix[center][i];
		}
		// the weight from center to itself is zero
		weights[center] = 0;
		for (int i = 0; i < num; i++) {
			// start from that center vertex, so skip when point to itself
			if (center == i) {
				continue;
			}
			int j = 0; // temporary index
			int k = 0; // the index of vertex which is minimum
			int min = INF; // find the minimum weight
			// find the minimum weight from the unsorted vertex
			while (j < num) {
				// if weights[j]=0, it mean j-th vertex has been added to the MST
				if (weights[j] != 0 && weights[j] < min) {
					min = weights[j];
					k = j;
				}
				j++;
			}
			// have the vertex that is minimum, add into the result array
			prims[index++] = mVexs[k];
			// label the weight=0, means the vertex has been added
			weights[k] = 0;
			// when the vertex is added to the result array, renew weight
			// weight of all vertex is based the the renewed center: k-th vertex
			for (j = 0; j < num; j++) {
				if (weights[j] != 0 && mMatrix[k][j] < weights[j]) {
					weights[j] = mMatrix[k][j];
				}
			}
		}
        // print MST
        for (int i = 0; i < index; i++) {
        	System.out.printf("%s ", prims[i]);
        }
        System.out.printf("\n");
	}

	private static class VERTEX {
		int weight; // 
		String end; // 
		
		public VERTEX(int weight, String end) {
			this.weight = weight;
			this.end = end;
		}
	}
	
	class MinHeap{
		private VERTEX [] data; 
	    private int size;    
	    private int maxSize; 

	    MinHeap(int maxSize){ // initialize
	        data = new VERTEX[maxSize+1];
	        data[0] = new VERTEX(Integer.MIN_VALUE, "V0");
	        this.size = 0;
	        this.maxSize = maxSize;
	    }
	    
	     public boolean insert(int weight, String inf){
	    	 VERTEX v = new VERTEX(weight, inf);
	         int index = ++size; 
	         if (size == maxSize){
	             System.out.println("堆已满，无法插入 "+weight);
	             return false;
	         }
	         while (weight < data[index/2].weight){
	             data[index] = data[index/2];
	             index = index/2;
	         }
	         data[index] = v;
	         return true;
	     }	
	}


	// MST using prim algorithm (binary heap)
	public void prim_BiHeap(int center) {
		int size = mMatrix.length;
		MinHeap heap = new MinHeap(size*size);
		int k = 0;
		int[] reWeights = new int[size];
		String[] reVerend = new String[size];
		VERTEX [] viewed = new VERTEX[size];
		int numView = 0;
		int j = 0;
		viewed[j++] = new VERTEX(0, label[center]);
		numView++;
		for(int n=0; n <size; n++) { // 
			for(int i = 0; i< size; i++) {
				heap.insert(mMatrix[center][i], label[i]);
			}
			VERTEX smallest = heap.data[0];
			int flag = 0;
			for(int p= 0; p < numView; p++) {
				if((viewed[p].end).equals(smallest.end)) {
					flag = 1; // has already viewed
					break;
				}
			}
			if(flag!=1) {
				viewed[j++] = smallest;
				numView++;
				reWeights[k++] = smallest.weight;
				reVerend[k++] = smallest.end;
				flag = 0;
			}
			char ch = smallest.end.charAt(1);
			center = (int)ch;
		}
		
		
	}

	// MST using Kruskal algorithm
	public void kruskal() {
		int index = 0; // index for array rets
		int[] vends = new int[mEdgNum]; // save end of every vertex
		EDGE[] rets = new EDGE[mEdgNum]; // result array, save edges of MST
		EDGE[] edges; // all edges of graph
		// get all the egdes from the graph
		edges = getEdges();
		// sort the edges according to the weights
		sortEdges(edges, mEdgNum);
		for (int i = 0; i < mEdgNum; i++) {
			int p1 = getPosition(edges[i].start); // get index of start of i-th edge
			int p2 = getPosition(edges[i].end); // get index of end of i-th edge
			int m = getEnd(vends, p1); // get end of p1
			int n = getEnd(vends, p2); // get end of p2
			if (m != n) { // edge i can be added, there is no circuit
				vends[m] = n;
				rets[index++] = edges[i]; // save result
			}
		}
        // print the information of MST
        int length = 0;
        for (int i = 0; i < index; i++) {
        	length += rets[i].weight;
        }
        System.out.printf("Kruskal=%d: ", length);
        for (int i = 0; i < index; i++) {
        	System.out.printf("(%s,%s) ", rets[i].start, rets[i].end);
        }
        System.out.printf("\n");
	}	
	
	// class for edges
	private static class EDGE {
		String start; // start of the edge
		String end; // end of the edge
		int weight; // weight of the edge

		public EDGE(String start, String end, int weight) {
			this.start = start;
			this.end = end;
			this.weight = weight;
		}
	}

	// get edges from the graph
	private EDGE[] getEdges() {
		int index = 0;
		EDGE[] edges;
		edges = new EDGE[mEdgNum];
		for (int i = 0; i < mVexs.length; i++) {
			for (int j = i + 1; j < mVexs.length; j++) {
				if (mMatrix[i][j] != INF) {
					edges[index++] = new EDGE(mVexs[i], mVexs[j], mMatrix[i][j]);
				}
			}
		}
		return edges;
	}

	// sort edges according to weights
	private void sortEdges(EDGE[] edges, int elen) {
		for (int i = 0; i < elen; i++) {
			for (int j = i + 1; j < elen; j++) {
				if (edges[i].weight > edges[j].weight) {
					// exchange edge i and edge j
					EDGE tmp = edges[i];
					edges[i] = edges[j];
					edges[j] = tmp;
				}
			}
		}
	}

	// get the end of i
	private int getEnd(int[] vends, int i) {
		while (vends[i] != 0) {
			i = vends[i];
		}
		return i;
	}

	static String[] label;

	// Random Connected Undirected Weighted Graph generator
	public static int[][] generator(int N, int S) { // N, size of the generated graph
		// S, sparseness (number of edges actually; from V-1 to V(V-1)/2)
		// random weight for the edge will range from 1 to 20
		// construct a adjacency matrix to store the graph
		int[][] RCUWG = new int[N][N];
		// generate label for graph
		String[] temp = new String[N];
		for (int i = 0; i < temp.length; i++) {
			temp[i] = "V" + (i + 1);
		}
		label = temp;
		// N nodes for sure, random the nodes sequence
		Integer[] index = new Integer[N];
		for (int i = 0; i < index.length; i++) {
			index[i] = i;
		}
		shuffle(index);
		// construct graph node by node, each time create a new node and a new edge
		Random df = new Random();
		int ranweigh = df.nextInt(20);
		for (int i = 1; i < index.length; i++) {
			// give the random weight
			ranweigh = df.nextInt(20);
			int tempInd = index[i];
			int tempPreInd = index[i - 1];
			RCUWG[tempPreInd][tempInd] = ranweigh + 1;
			RCUWG[tempInd][tempPreInd] = ranweigh + 1;
		}
		// if there are S-N+1 nodes left, randomly give edges
		int leftNode = S - N + 1;
		for (int i = 1; i < leftNode + 1; i++) {
			int connect1 = df.nextInt(N - 1); // choose a random node between 0 and N-1
			int connect2 = df.nextInt(N - 1); // choose another
			// if meet two same nodes by accident, then continue one more time
			if (connect1 == connect2 || RCUWG[connect1][connect2] != 0) {
				i--;
				continue;
			}
			// give the random weight
			ranweigh = df.nextInt(20);
			RCUWG[connect1][connect2] = ranweigh + 1;
			RCUWG[connect2][connect1] = ranweigh + 1;
		}
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				if (RCUWG[i][j] == 0 && i != j) {
					RCUWG[i][j] = Integer.MAX_VALUE;
				}
			}
		}
		return RCUWG;
	}

	// for randomly shuffle the array number, which means randomly connect with node
	// one by one
	private static Random rand = new Random();

	public static <T> void swap(T[] a, int i, int j) {
		T temp = a[i];
		a[i] = a[j];
		a[j] = temp;
	}

	public static <T> void shuffle(T[] arr) {
		int length = arr.length;
		for (int i = length; i > 0; i--) {
			int randInd = rand.nextInt(i);
			swap(arr, randInd, i - 1);
		}
	}

	public static void main(String[] args) {
		int[][] graph = generator(10, 30);
		MST_Algorithm pG = new MST_Algorithm(label, graph);
        pG.print(); // print graph in standard
        pG.prim_Array(0);
        pG.kruskal();
		// compare the efficiency
//		for (int i = 10; i < 3000; i++) {
//			//graph = generator(i, (int)(i*(i-1)*0.8/2)); // dense graph creation
//			graph = generator(i, 2*i); // sparse graph creation
//			pG = new MST_Algorithm(label, graph);
//			long startTime = System.currentTimeMillis();
//			 pG.kruskal(); // kruskal
//			// pG.prim_BiHeap(0);
//			//pG.prim_Array(0); // prim algorithm using array
//			long endTime = System.currentTimeMillis();
//			long consume = endTime - startTime;
//			System.out.println(consume);
//		}

	}
}


