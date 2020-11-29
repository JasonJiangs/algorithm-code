package as12;

import java.util.ArrayList;
import java.util.NoSuchElementException;

public class BellmanFordVsDijkstra {

	public static void main(String[] args) {
		long time1;
		long time2;
		long time3;
		long time4;
		long start1;
		long end1;
		long start2;
		long end2;
		long start3;
		long end3;
		long start4;
		long end4;
		int start = 100;
		int end = 5000;
		int gap = 100;
		int rep = 1;
		Graph g;

//		g=new Graph(10);
//		g.generateDense();
////		g.generateSparse();
//		g.print();
//		BellmanFord bf=new BellmanFord(g);
//		bf.solve(2);
//		Dijkstra dj=new Dijkstra(g);
//		dj.solveArray(2);
//		dj.solveBinaryHeap(2);
//		dj.solveFibHeap(2);

		for (int i = start; i < end; i += gap) {
			time1 = 0;// BellmanFord
			time2 = 0;// Dijkstra with array
			time3 = 0;// Dijkstra with binaryHeap
			time4 = 0;// Dijkstra with FibnacciHeap
			g = new Graph(i);
//			g.generateDense();
			g.generateSparse();
			BellmanFord bf = new BellmanFord(g);
			Dijkstra dj = new Dijkstra(g);
			for (int j = 0; j < rep; j++) {
				start1 = System.currentTimeMillis();
				bf.solve(0);
				end1 = System.currentTimeMillis();
//				start2 = System.currentTimeMillis();
//				dj.solveArray(0);
//				end2 = System.currentTimeMillis();
//				start3 = System.currentTimeMillis();
//				dj.solveBinaryHeap(0);
//				end3 = System.currentTimeMillis();
//				start4 = System.currentTimeMillis();
//				dj.solveFibHeap(0);
//				end4 = System.currentTimeMillis();
				time1 += (end1 - start1);
				//time2 += (end2 - start2);
				//time3 += (end3 - start3);
//				time4 += (end4 - start4);
			}
			System.out.println(time1 / rep);
		}
//		 "       " + time1 / rep + "        " + time2 / rep + "         " + time3 / rep
//			+ "         " + 
	}

}

class Dijkstra {
	Graph g;
	boolean s[];// if in the result
	int[] d;
	BinaryHeap bh;
	FibHeap fh;

	Dijkstra(Graph g) {
		this.g = g;
		s = new boolean[g.nbnodes];
		d = new int[g.nbnodes];
		bh = new BinaryHeap(g.nbnodes);
		fh = new FibHeap();
	}

	public void relax(Node u, Node v, int weight) {
		if (!s[v.index] && d[u.index] != Integer.MAX_VALUE && d[v.index] > d[u.index] + weight) {
			d[v.index] = d[u.index] + weight;
		}
	}

	public void relaxInBinaryHeap(Node u, Node v, int weight) {
		if (!v.ifInResult && u.key != Integer.MAX_VALUE && v.key > u.key + weight) {
			for (int j = 0; j < bh.heapSize; j++) {
				if (bh.heap[j] == v) {
					bh.delete(j);
					v.key = u.key + weight;
					bh.insert(v);
					break;
				}
			}
		}
	}

	public void relaxInFibHeap(Node u, Node v, int weight) {
		if (!v.ifInResult && u.key != Integer.MAX_VALUE && v.key > u.key + weight) {
			fh.decreaseKey(v, u.key + weight);
		}
	}

	public int minDistanceIndex(int[] d, boolean[] s) {// choose the minimum in an array
		int min = Integer.MAX_VALUE;
		int minIndex = -1;
		for (int i = 0; i < g.nbnodes; i++)
			if (s[i] == false && d[i] <= min) { // if i in the array and d[i]<min
				min = d[i];
				minIndex = i;
			}
		return minIndex;
	}

	public void solveArray(int source) {
		for (int i = 0; i < g.nbnodes; i++) {// Initialization
			d[i] = Integer.MAX_VALUE;
			s[i] = false;
		}
		d[source] = 0;// root
		for (int i = 0; i < g.nbnodes; i++) {// v times
			int u = minDistanceIndex(d, s);// extract-min in the array
			s[u] = true;// put the min in the result, remove from the array
			Node nu = g.nodes.get(u);

			for (Edge e : nu.adjecentEdge) {// for all the minimum one's adjecentEdge
				Node nv = e.destination;
				int weight = e.weight;
				relax(nu, nv, weight);// relax adjecentEdge
			}
		}

//		System.out.println("Dijkstra With Array: Vertex Distance from Source " + source);
//		for (int i = 0; i < g.nbnodes; i++) {
//			System.out.println(source + "--->" + i + " = " + d[i]);
//		}

	}

	public void solveBinaryHeap(int source) {
		g.nodes.get(source).key = 0;// root
		g.nodes.get(source).ifInResult = false;
		bh.insert(g.nodes.get(source));
		for (int i = 0; i < g.nbnodes; i++) {// Initialization
			if (i != source) {
				g.nodes.get(i).key = Integer.MAX_VALUE;// d[i]=Integer.MAX_VALUE
				g.nodes.get(i).ifInResult = false;// s[i]=false
				bh.insert(g.nodes.get(i));
			}
		}
		for (int i = 0; i < g.nbnodes; i++) {
			Node u = bh.deleteMin();// extract-min in the heap and remove it from the heap
			u.ifInResult = true;

			for (Edge e : u.adjecentEdge) {// for all the minimum one's adjecentEdge
				Node v = e.destination;
				int weight = e.weight;
				relaxInBinaryHeap(u, v, weight);
			}
		}

//		System.out.println("Dijkstra With BinaryHeap: Vertex Distance from Source " + source);
//		for (int i = 0; i < g.nbnodes; i++) {
//			System.out.println(source + "--->" + i + " = " + g.nodes.get(i).key);
//		}

	}

	public void solveFibHeap(int source) {
		g.nodes.get(source).ifInResult = false;
		fh.insert(g.nodes.get(source), 0);// root
		for (int i = 0; i < g.nbnodes; i++) {// Initialization
			if (i != source) {
				g.nodes.get(i).ifInResult = false;// s[i]=false
				fh.insert(g.nodes.get(i), Integer.MAX_VALUE); // d[i]=Integer.MAX_VALUE
			}
		}
		for (int i = 0; i < g.nbnodes; i++) {
			Node u = fh.removeMin();// extract-min in the heap and remove it from the heap
			u.ifInResult = true;
			for (Edge e : u.adjecentEdge) {// for all the minimum one's adjecentEdge
				Node v = e.destination;
				int weight = e.weight;
				relaxInFibHeap(u, v, weight);
			}
		}

//		System.out.println("Dijkstra With FibHeap: Vertex Distance from Source " + source);
//		for (int i = 0; i < g.nbnodes; i++) {
//			System.out.println(source + "--->" + i + " = " + g.nodes.get(i).key);
//		}

	}

}

class BellmanFord {
	Graph g;
	int[] d;
	Node[] p;// predecessor

	BellmanFord(Graph g) {
		this.g = g;
		d = new int[g.nbnodes];
		p = new Node[g.nbnodes];
	}

	public void relax(Node u, Node v, int weight) {
		if (d[u.index] != Integer.MAX_VALUE && d[v.index] > d[u.index] + weight) {
			d[v.index] = d[u.index] + weight;
			p[v.index] = u;
		}
	}

	public boolean detectNegativeCycles(Edge e) {// if d[v]>d[u]+weight have cycle
		if (d[e.destination.index] > d[e.source.index] + e.weight) {
			return false;
		}
		return true;
	}

	public void solve(int source) {
		for (int i = 0; i < g.nbnodes; i++) {// 1.Initialize-Single-Source(V,s)
			d[i] = Integer.MAX_VALUE;
			p[i] = null;
		}
		d[source] = 0;// start from the root
		for (int i = 0; i < g.nbnodes - 1; i++) {// each edge relax v-1 times
			for (int j = 0; j < g.nbEdge; j++) {
				relax(g.edges.get(j).source, g.edges.get(j).destination, g.edges.get(j).weight);// Relax(u,v,w)
			}
		}
		for (int i = 0; i < g.nbEdge; i++) {// detect cycle
			if (!detectNegativeCycles(g.edges.get(i)) && d[g.edges.get(i).source.index] != Integer.MAX_VALUE) {
				System.out.println("Graph contains negative weight cycle");
			}
		}

//		System.out.println("BellmanFord: Bertex Distance from Source " + source);
//		for (int i = 0; i < g.nbnodes; i++) {
//			System.out.println(source + "--->" + i + " = " + d[i]);
//		}

	}
}

class Graph {
	public int nbnodes;
	public int nbEdge;
	public ArrayList<Node> nodes;
	public ArrayList<Edge> edges;

	public Graph(int nbnodes) {
		this.nbnodes = nbnodes;
		this.nodes = new ArrayList<>();
		this.edges = new ArrayList<>();
	}

	public void generateSparse() {
		nbEdge = 0;
		Node first = new Node(0);
		nodes.add(first);
		Node current;
		Node previous;
		for (int i = 1; i < nbnodes; i++) {
			current = new Node(i);
			previous = nodes.get(i - 1);
			nodes.add(current);
			addEdge(current, previous);
			nbEdge += 2;
		}
	}

	public void generateDense() {
		nbEdge = 0;
		Node first = new Node(0);
		nodes.add(first);
		Node previous;
		for (int i = 1; i < nbnodes; i++) {
			Node current = new Node(i);
			nodes.add(current);
			for (int j = 0; j < i; j++) {
				previous = nodes.get(j);
				addEdge(current, previous);
				nbEdge += 2;
			}
		}
	}

	public void print() {
		System.out.println("|V|=" + nbnodes + " |E|=" + nbEdge);
		for (int i = 0; i < nbnodes; i++) {
			System.out.print(i + " --> {   ");
			for (int j = 0; j < nodes.get(i).adjecentEdge.size(); j++) {
				System.out.print("[ " + nodes.get(i).adjecentEdge.get(j).destination.index + ":"
						+ nodes.get(i).adjecentEdge.get(j).weight + " ] ,  ");
			}
			System.out.println("}");
		}
	}

	public void addEdge(Node current, Node previous) {
		Edge e;
		int r1 = (int) (Math.random() * 101);
		int r2 = (int) (Math.random() * 101);
		e = new Edge(current, previous, r1);
		current.adjecentEdge.add(e);
		edges.add(e);
		e = new Edge(previous, current, r2);
		previous.adjecentEdge.add(e);
		edges.add(e);
	}
}

class Node {

	public int index;
	public int key;// distance
	ArrayList<Edge> adjecentEdge;
	public boolean ifInResult;
	int degree;
	Node left;
	Node right;
	Node child;
	Node parent;
	boolean mark;

	public Node(int index) {
		this.index = index;
		adjecentEdge = new ArrayList<>();
		this.degree = 0;
		this.mark = false;
		this.left = this;
		this.right = this;
		this.parent = null;
		this.child = null;
	}
}

class Edge implements Comparable<Edge> {
	Node source;
	Node destination;
	int weight;

	public Edge(Node source, Node destination, int weight) {
		this.source = source;
		this.destination = destination;
		this.weight = weight;
	}

	public int compareTo(Edge compareEdge) {
		return this.weight - compareEdge.weight;
	}
}

class BinaryHeap {
	/** The number of children each node has **/
	private static final int d = 2;
	public int heapSize;
	public Node[] heap;

	/** Constructor **/
	public BinaryHeap(int capacity) {
		heapSize = 0;
		heap = new Node[capacity + 1];
	}

	/** Function to check if heap is empty **/
	public boolean isEmpty() {
		return heapSize == 0;
	}

	/** Check if heap is full **/
	public boolean isFull() {
		return heapSize == heap.length;
	}

	/** Clear heap */
	public void makeEmpty() {
		heapSize = 0;
	}

	/** Function to get index parent of i **/
	private int parent(int i) {
		return (i - 1) / d;
	}

	/** Function to get index of k th child of i **/
	private int kthChild(int i, int k) {
		return d * i + k;
	}

	/** Function to insert element */
	public void insert(Node x) {
		if (isFull())
			throw new NoSuchElementException("Overflow Exception");
		/** Percolate up **/
		heap[heapSize++] = x;
		heapifyUp(heapSize - 1);
	}

	/** Function to find least element **/
	public Node findMin() {
		if (isEmpty())
			throw new NoSuchElementException("Underflow Exception");
		return heap[0];
	}

	/** Function to delete min element **/
	public Node deleteMin() {
		Node keyItem = heap[0];
		delete(0);
		return keyItem;
	}

	/** Function to delete element at an index **/
	public Node delete(int ind) {
		if (isEmpty())
			throw new NoSuchElementException("Underflow Exception");
		Node keyItem = heap[ind];
		heap[ind] = heap[heapSize - 1];
		heapSize--;
		heapifyDown(ind);
		return keyItem;
	}

	/** Function heapifyUp **/
	private void heapifyUp(int childInd) {
		Node tmp = heap[childInd];
		while (childInd > 0 && tmp.key < heap[parent(childInd)].key) {
			heap[childInd] = heap[parent(childInd)];
			childInd = parent(childInd);
		}
		heap[childInd] = tmp;
	}

	/** Function heapifyDown **/
	public void heapifyDown(int ind) {
		int child;
		Node tmp = heap[ind];
		while (kthChild(ind, 1) < heapSize) {
			child = minChild(ind);
			if (heap[child].key < tmp.key)
				heap[ind] = heap[child];
			else
				break;
			ind = child;
		}
		heap[ind] = tmp;
	}

	/** Function to get smallest child **/
	private int minChild(int ind) {
		int bestChild = kthChild(ind, 1);
		int k = 2;
		int pos = kthChild(ind, k);
		while ((k <= d) && (pos < heapSize)) {
			if (heap[pos].key < heap[bestChild].key)
				bestChild = pos;
			pos = kthChild(ind, k++);
		}
		return bestChild;
	}

	/** Function to print heap **/
	public void printHeap() {
		System.out.print("\nHeap = ");
		for (int i = 0; i < heapSize; i++)
			System.out.print(heap[i].index + "." + heap[i].key + "  ");
		System.out.println();
	}
}

class FibHeap {

	private Node minNode;
	private int nnodes;

	public FibHeap() // Fibonacci heap default constructor
	{

	}

	public boolean isEmpty() // Returns null if heap is empty
	{
		return minNode == null;
	}

	public void clear() // Resets the heap
	{
		minNode = null;
		nnodes = 0;
	}

	// Decalration to keep the fibonacci heap property
	private static final double oneOverLogPhi = 1.0 / Math.log((1.0 + Math.sqrt(5.0)) / 2.0);

	// insert functions inserts the node in the heap in O(1) time.
	public Node insert(Node node, int key)

	{
		node.key = key;
		if (minNode != null) {
			node.left = minNode;
			node.right = minNode.right;
			minNode.right = node;
			node.right.left = node;

			if (key < minNode.key) {
				minNode = node;
			}

		} else {
			minNode = node;
		}

		nnodes++;
		return node;
	}

	public Node min() // This function returns Minimum node in the heap.
	{
		return minNode;
	}

	public int size() // This function returns the current size of the haep
	{
		return nnodes;
	}

	/*
	 * This function is a helper function for consolidate function for linking of
	 * node y and node x in the heap
	 */
	protected void link(Node y, Node x) {
		// remove y from root list of heap
		y.left.right = y.right;
		y.right.left = y.left;

		// make y a child of x
		y.parent = x;

		if (x.child == null) {
			x.child = y;
			y.right = y;
			y.left = y;
		} else {
			y.left = x.child;
			y.right = x.child.right;
			x.child.right = y;
			y.right.left = y;
		}

		// increase degree[x]
		x.degree++;

		// set mark[y] false
		y.mark = false;
	}

	/**
	 * Removes the smallest element from the heap. This will cause the trees in the
	 * heap to be consolidated, if necessary. Running time is O(log n)
	 *
	 * This function returns node with the smallest key
	 */
	public Node removeMin() {
		Node z = minNode;

		if (z != null) {
			int numKids = z.degree;
			Node x = z.child;
			Node tempRight;

			// for each child of z do...
			while (numKids > 0) {
				tempRight = x.right;

				// remove x from child list
				x.left.right = x.right;
				x.right.left = x.left;

				// add x to root list of heap
				x.left = minNode;
				x.right = minNode.right;
				minNode.right = x;
				x.right.left = x;

				// set parent[x] to null
				x.parent = null;
				x = tempRight;
				numKids--;
			}

			// remove z from root list of heap
			z.left.right = z.right;
			z.right.left = z.left;

			if (z == z.right) {
				minNode = null;
			} else {
				minNode = z.right;
				consolidate();
			}

			// decrement size of heap
			nnodes--;
		}

		return z;
	}

	protected void consolidate() {
		int arraySize = ((int) Math.floor(Math.log(nnodes) * oneOverLogPhi)) + 1;

		ArrayList<Node> array = new ArrayList<Node>(arraySize);

		// Initialize degree array
		for (int i = 0; i < arraySize; i++) {
			array.add(null);
		}

		// Find the number of root nodes.
		int numRoots = 0;
		Node x = minNode;

		if (x != null) {
			numRoots++;
			x = x.right;

			while (x != minNode) {
				numRoots++;
				x = x.right;
			}
		}

		// For each node in root list do...
		while (numRoots > 0) {
			// Access this node's degree..
			int d = x.degree;
			Node next = x.right;

			// ..and see if there's another of the same degree.
			for (;;) {
				Node y = array.get(d);
				if (y == null) {
					// Nope.
					break;
				}

				// There is, make one of the nodes a child of the other.
				// Do this based on the key value.
				if (x.key > y.key) {
					Node temp = y;
					y = x;
					x = temp;
				}

				// Node<T> y disappears from root list.
				link(y, x);

				// We've handled this degree, go to next one.
				array.set(d, null);
				d++;
			}

			// Save this node for later when we might encounter another
			// of the same degree.
			array.set(d, x);

			// Move forward through list.
			x = next;
			numRoots--;
		}

		// Set min to null (effectively losing the root list) and
		// reconstruct the root list from the array entries in array[].
		minNode = null;

		for (int i = 0; i < arraySize; i++) {
			Node y = array.get(i);
			if (y == null) {
				continue;
			}

			// We've got a live one, add it to root list.
			if (minNode != null) {
				// First remove node from root list.
				y.left.right = y.right;
				y.right.left = y.left;

				// Now add to root list, again.
				y.left = minNode;
				y.right = minNode.right;
				minNode.right = y;
				y.right.left = y;

				// Check if this is a new min.
				if (y.key < minNode.key) {
					minNode = y;
				}
			} else {
				minNode = y;
			}
		}
	}

	/**
	 * Performs a cascading cut operation. This cuts y from its parent and then does
	 * the same for its parent, and so on up the tree.
	 *
	 * Running time is: O(log n);
	 */

	protected void cascadingCut(Node y) {
		Node z = y.parent;

		// if there's a parent...
		if (z != null) {
			// if y is unmarked, set it marked
			if (!y.mark) {
				y.mark = true;
			} else {
				// it's marked, cut it from parent
				cut(y, z);

				// cut its parent as well
				cascadingCut(z);
			}
		}
	}

	/**
	 * Decreases the key value for a heap node, given the new value to take on. The
	 * structure of the heap may be changed and will not be consolidated.
	 *
	 * Running time: O(1) amortized
	 */
	public void decreaseKey(Node x, int k) {
		if (k > x.key) {
			return;
			// throw new IllegalArgumentException(
			// "decreaseKey() got larger key value");
		}

		x.key = k;

		Node y = x.parent;

		if ((y != null) && (x.key < y.key)) {
			cut(x, y);
			cascadingCut(y);
		}

		if (x.key < minNode.key) {
			minNode = x;
		}
	}

	/*
	 * The reverse of the link operation: removes x from the child list of y. This
	 * method assumes that min is non-null.
	 *
	 * Running time is: O(1)
	 */

	protected void cut(Node x, Node y) {
		// remove x from childlist of y and decrement degree[y]
		x.left.right = x.right;
		x.right.left = x.left;
		y.degree--;

		// reset y.child if necessary
		if (y.child == x) {
			y.child = x.right;
		}

		if (y.degree == 0) {
			y.child = null;
		}

		// add x to root list of heap
		x.left = minNode;
		x.right = minNode.right;
		minNode.right = x;
		x.right.left = x;

		// set parent[x] to nil
		x.parent = null;

		// set mark[x] to false
		x.mark = false;
	}
}
