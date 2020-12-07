package as13;

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Deque;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.NoSuchElementException;
import java.util.Scanner;
import java.util.Set;

public class FordFulkerson_Maxflow {

	public static void main(String[] args) {
		Scanner input = new Scanner(System.in);
		FlowNetwork<Integer> fn = new FlowNetwork<Integer>();
		System.out.println("Vertex number:");
		int vnum = input.nextInt();
		// for dense garph
//		linkedNode[] list = generator.listGenerator(vnum, (int)((vnum))*(vnum-1));//|V|*|V-1|
		// for spase graph
		linkedNode[] list = generator.listGenerator(vnum, (int) ((vnum - 1))); // |V|
		int[][] matrix = new int[vnum][vnum];
		matrix = convertor.listToMatrix(list);
		IntegralDirectedGraph<Integer> graph = new IntegralDirectedGraph<Integer>();
		graph.addNode(1);
		for (int i = 0; i < vnum; i++) {
			if (i == 0) {
				for (int j = 0; j < vnum; j++) {
					graph.addNode(j);
				}
			}
			for (int j = 0; j < vnum; j++) {
				graph.addEdge(i, j, matrix[i][j]);
				System.out.print(matrix[i][j] + " ");
			}
			System.out.println();
		}

		// DFS time consuming and maxflow result
		long timedfs, timebfs;
		timedfs = System.nanoTime();
		fn = Fordfulkenson_Generic.maxFlow(graph, 0, (vnum - 1));
		Fordfulkenson_Generic.findMaxFlow(fn, 0, (vnum - 1));
		int maxflow = 0;
		for (int i = 1; i < vnum; i++) {
			maxflow += fn.getEdge(0, i).getFlow();
		}
		timedfs = System.nanoTime() - timedfs;
		System.out.println("DFS version's maxflow: " + maxflow + " & Time consuming: " + timedfs);

		// BFS time consuming and maxflow result
		timebfs = System.nanoTime();
		fn = FordFulkersonScaling.maxFlow(graph, 0, (vnum - 1));
		FordFulkersonScaling.findMaxFlow(fn, 0, (vnum - 1));
		maxflow = 0;
		for (int i = 1; i < vnum; i++) {
			maxflow += fn.getEdge(0, i).getFlow();
		}
		timebfs = System.nanoTime() - timebfs;
		System.out.println("BFS version's maxflow: " + maxflow + " & Time consuming: " + timebfs);
	}
}


class generator {
	public static final int maxWeight = 20; // 100
	private static boolean inputValidity(int size, int sparseness) {
		if ((sparseness < size - 1) || (sparseness > (size - 1) * size) || size < 0) {
			return false;
		} else {
			return true;
		}
	}
	// Random Connected Directed Weighted Graph Generator --> AdjacentMatrix
	public static int[][] matrixGenerator(int size, int sparseness) {
		if (!inputValidity(size, sparseness))
			return null;
		// Create the adjacentMatrix
		int[][] matrix = new int[size][size]; // --> All the item is initialized to 0 in default

		// Make the graph connected
		int anotherPoint, weight; // Another Point
		for (int index = 1; index < size; index++) {
			// Generate another random point
			anotherPoint = (int) (Math.random() * index);
			weight = 1 + (int) (Math.random() * size);
			weight = 1 + (int) (Math.random() * maxWeight);

			// Link the 2 points in Small --> Large Direction
			matrix[anotherPoint][index] = weight;
		}
		sparseness = sparseness - size + 1;
		// Deal with the sparseness: Create random connection P1 --> P2
		int point_1, point_2;
		while (sparseness > 0) {
			// Generate 2 different point's Index
			point_1 = (int) (Math.random() * size);
			do {
				point_2 = (int) (Math.random() * size);
			} while (point_1 == point_2);
			// Insert the edge
			if (matrix[point_1][point_2] == 0) {
				// Insert if no edge exists
				matrix[point_1][point_2] = (int) (Math.random() * size);
				sparseness--;
			}
		}
		return matrix;
	}

	public static void printMatrix(int[][] matrix) {
		int size = matrix.length;
		// Print the first line: NodeIndex
		System.out.printf("%4s", "");
		for (int index = 0; index < size; index++) {
			System.out.printf("%4d", index);
		}
		System.out.print("\n");
		// Print the following line
		for (int index = 0; index < size; index++) {
			// First is NodeIndex
			System.out.printf("%4d", index);
			for (int item : matrix[index]) {
				System.out.printf("%4d", item);
			}
			System.out.print("\n");
		}
	}

	// Random Connected Directed Weighted Graph Generator --> AdjacentList
	public static linkedNode[] listGenerator(int size, int sparseness) {
		// Check the input Validity
		if (!inputValidity(size, sparseness))
			return null;
		// Create the adjacentList
		linkedNode[] list = new linkedNode[size];
		// Initialize the adjacentList
		for (int index = 0; index < size; index++) {
			list[index] = new linkedNode(index, 0);
		}
		// Make the graph connected
		int anotherPoint, weight; // Another Point
		for (int index = 1; index < size; index++) {
			// Generate another random point
			anotherPoint = (int) (Math.random() * index);
			// weight = 1 + (int) (Math.random() * size);
			weight = 1 + (int) (Math.random() * maxWeight);

			// AnotherPoint --> Index
			list[anotherPoint].tail().next = new linkedNode(index, weight);
		}
		sparseness = sparseness - size + 1;
		// Deal with the sparseness: Create random connection P1 --> P2
		int point_1, point_2;
		while (sparseness > 0) {
			// Generate 2 different point's Index
			point_1 = (int) (Math.random() * size);
			do {
				point_2 = (int) (Math.random() * size);
			} while (point_1 == point_2);
			// Check whether the edge exists
			linkedNode pointer = list[point_1];
			boolean exist = false;
			while (pointer.next != null) {
				if (pointer.next.id == point_2) {
					exist = true;
					break;
				}
				pointer = pointer.next;
			}
			// Insert the edge if the connection not Exist
			if (!exist) {
				pointer.next = new linkedNode(point_2, 1 + (int) (Math.random() * maxWeight));
				sparseness--;
			} else {
				continue;
			}
		}
		return list;
	}

	public static void printList(linkedNode[] list) {
		for (linkedNode root : list) {
			// Print the Root Node
			System.out.print(root + " -> ");
			root = root.next;
			while (root != null) {
				System.out.printf("%7s", root);
				root = root.next;
			}
			System.out.print("\n");
		}
	}

	public static void main(String[] args) {
		printList(listGenerator(5, 9));
	}
}

//A linkedList used for adjacentList
class linkedNode {
	@Override
	public String toString() {
		return this.id + "(" + this.weight + ")";
	}
	linkedNode next;
	int weight;
	int id;
	linkedNode(int id, int weight) {
		this.id = id;
		this.weight = weight;
		this.next = null;
	}
	linkedNode(int id, int weight, linkedNode next) {
		this.id = id;
		this.weight = weight;
		this.next = next;
	}
	linkedNode tail() {
		linkedNode pointer = this;
		while (!(pointer.next == null)) {
			pointer = pointer.next;
		}
		return pointer;
	}
}

final class IntegralDirectedGraph<T> implements Iterable<T> {
//	 A map from nodes in the graph to sets of outgoing edges. Each set of edges is
//	 represented by a map from edges to doubles.
	private final Map<T, Map<T, Integer>> mGraph = new HashMap<T, Map<T, Integer>>();
	public boolean addNode(T node) {
		//If the node already exists, don't do anything. 
		if (mGraph.containsKey(node))
			return false;
		//Otherwise, add the node with an empty set of outgoing edges.
		mGraph.put(node, new HashMap<T, Integer>());
		return true;
	}

	public void addEdge(T start, T dest, int length) {
		// Confirm both endpoints exist.
		if (!mGraph.containsKey(start) || !mGraph.containsKey(dest))
			throw new NoSuchElementException("Both nodes must be in the graph.");
		// Add the edge. 
		mGraph.get(start).put(dest, length);
	}
	
	public void removeEdge(T start, T dest) {
		// Confirm both endpoints exist.
		if (!mGraph.containsKey(start) || !mGraph.containsKey(dest))
			throw new NoSuchElementException("Both nodes must be in the graph.");
		mGraph.get(start).remove(dest);
	}

	public Map<T, Integer> edgesFrom(T node) {
		// Check that the node exists.
		Map<T, Integer> arcs = mGraph.get(node);
		if (arcs == null)
			throw new NoSuchElementException("Source node does not exist.");
		return Collections.unmodifiableMap(arcs);
	}

	public Iterator<T> iterator() {
		return mGraph.keySet().iterator();
	}

	public boolean containsNode(T node) {
		return mGraph.containsKey(node);
	}

	public int size() {
		return mGraph.size();
	}

	public boolean isEmpty() {
		return mGraph.isEmpty();
	}
}


class convertor {
    // Random Connected Directed Weighted Graph Convertor --> Convert the AdjacentMatrix
    public static linkedNode[] matrixToList(int[][] generatedMatrix) {
        int size = generatedMatrix.length;
        linkedNode[] list = new linkedNode[size];
        // Initialize the adjacentList
        for (int index = 0; index < size; index++) {
            list[index] = new linkedNode(index, 0);
        }
        // Access each item and Do insertion
        linkedNode pointer;
        int arrayIndex = 0, itemIndex;
        for (int[] array: generatedMatrix) {
            itemIndex = 0;
            pointer = list[arrayIndex++];
            for(int item: array){
                if(item != 0){
                    pointer.next = new linkedNode(itemIndex, item);
                    pointer = pointer.next;
                }
                itemIndex++;
            }
        }
        return list;
    }

    public static int[][] listToMatrix(linkedNode[] list){
        // Form the Result Matrix
        int length = list.length;
        int[][] matrix = new int[length][length];
        // Loop all the items in the List
        int rootIndex = 0;
        for(linkedNode pointer: list){
            pointer = pointer.next;
            while(pointer != null){
                matrix[rootIndex][pointer.id] = pointer.weight;
                pointer = pointer.next;
            }
            rootIndex = rootIndex + 1;
        }
        return matrix;
    }
}


final class Fordfulkenson_Generic {
    public static <T> FlowNetwork<T> maxFlow(IntegralDirectedGraph<T> g, T s, T t) {
        // Construct the structure of the resulting flow network.
        FlowNetwork<T> result = new FlowNetwork<T>();
        // Copy over nodes. 
        for (T node: g)
            result.addNode(node);
        // Copy over edges.
        for (T node: g)
            for (Map.Entry<T, Integer> edge: g.edgesFrom(node).entrySet())
                result.addEdge(node, edge.getKey()).setCapacity(edge.getValue());
        /* Compute a max-flow in this flow network. */
        findMaxFlow(result, s, t);
        return result;
    }

    public static <T> void findMaxFlow(FlowNetwork<T> g, T s, T t) {
        /* Confirm that s and t are valid. */
        if (!g.containsNode(s) || !g.containsNode(t))
            throw new NoSuchElementException("Start and end nodes must be in hthe flow network!");
        if (isEqual(s, t)) return;
        ResidualGraph<T> gResidual = new ResidualGraph<T>(g);
        while (true) {
            /* Find an augmenting s-t path in the residual graph. */
            Deque<ResidualGraph.Edge<T>> path = findPath(s, t, gResidual);
            if (path == null) break;
            augmentPath(path);
        }
        for (T node: gResidual) {
            for (ResidualGraph.Edge<T> edge: gResidual.edgesFrom(node)) {
                if (!edge.isOriginal())
                    g.getEdge(edge.getEnd(), edge.getStart()).setFlow(edge.getCapacity());
            }
        }
    }

    private static <T> boolean isEqual(T one, T two) {
        /* If either are null, they're equal only if they're both null. */
        if (one == null || two == null)
            return one == null && two == null;
        return one.equals(two);
    }

    private static <T> Deque<ResidualGraph.Edge<T>> findPath(T start, T dest,
                                              ResidualGraph<T> graph) {
        return findPathRec(start, dest, graph, new HashSet<T>());
    }
    
    private static <T> Deque<ResidualGraph.Edge<T>> findPathRec(T start, T dest,
                                                 ResidualGraph<T> graph,
                                                 Set<T> visited) {
        if (visited.contains(start)) return null;
        visited.add(start);
        if (isEqual(start, dest))
            return new ArrayDeque<ResidualGraph.Edge<T>>();
        for (ResidualGraph.Edge<T> edge: graph.edgesFrom(start)) {
            if (edge.getCapacity() == 0) continue;
            Deque<ResidualGraph.Edge<T>> result = findPathRec(edge.getEnd(), dest,
                                               graph, visited);
            if (result != null) {
                result.addFirst(edge);
                return result;
            }
        }
        return null;
    }
    private static <T> void augmentPath(Deque<ResidualGraph.Edge<T>> path) {
        int capacity = Integer.MAX_VALUE;
        for (ResidualGraph.Edge<T> edge: path)
            capacity = Math.min(capacity, edge.getCapacity());
        for (ResidualGraph.Edge<T> edge: path)
            edge.addFlow(capacity);
    }
}


final class FordFulkersonScaling {
    public static <T> FlowNetwork<T> maxFlow(IntegralDirectedGraph<T> g, T s, T t) {
        /* Construct the structure of the resulting flow network. */
        FlowNetwork<T> result = new FlowNetwork<T>();
        /* Copy over nodes. */
        for (T node: g)
            result.addNode(node);
        /* Copy over edges. */
        for (T node: g)
            for (Map.Entry<T, Integer> edge: g.edgesFrom(node).entrySet())
                result.addEdge(node, edge.getKey()).setCapacity(edge.getValue());
        /* Compute a max-flow in this flow network. */
        findMaxFlow(result, s, t);
        return result;
    }

    public static <T> void findMaxFlow(FlowNetwork<T> g, T s, T t) {
        Map<FlowNetwork.Edge<T>, Integer> edges = new HashMap<FlowNetwork.Edge<T>, Integer>();
        for (T node: g) {
            for (FlowNetwork.Edge<T> edge: g.edgesFrom(node)) {
                /* Cache the capacity for later on. */
                edges.put(edge, edge.getCapacity());
                edge.setFlow(0);
                edge.setCapacity(0);
            }
        }

        /* Run the capacity-scaling rounds. */
        for (int bit = numRequiredBits(edges.values()); bit >= 0; --bit) {
            /* Scan across all edges, uncovering the next bit. */
            for (Map.Entry<FlowNetwork.Edge<T>, Integer> edge: edges.entrySet()) {
                edge.getKey().setCapacity(2 * edge.getKey().getCapacity() +
                                          ((edge.getValue() & (1 << bit)) >>> bit));
                edge.getKey().setFlow(2 * edge.getKey().getFlow());
            }
            /* Run another iteration of Ford-Fulkerson on this flow graph. */
            Fordfulkenson_Generic.findMaxFlow(g, s, t);
        }
    }

    private static int numRequiredBits(Collection<Integer> capacities) {
        int numBits = -1;
        for (Integer capacity: capacities)
            numBits = Math.max(numBits, 32 - Integer.numberOfLeadingZeros(capacity));
        return numBits;
    }
}


final class FlowNetwork<T> implements Iterable<T> {
    public static final class Edge<T> {
        private final T start; // Node start point
        private final T end;   // Node endpoint
        private int capacity;  // The capacity of the edge
        private int flow;      // The flow across the edge.

        public T getStart() {
            return start;
        }

        public T getEnd() {
            return end;
        }

        public int getCapacity() {
            return capacity;
        }

        public void setCapacity(int capacity) {
            /* Check that this is a valid capacity. */
            if (capacity < 0)
                throw new IllegalArgumentException("Capacities must be non-negative.");
            if (flow > capacity)
                throw new IllegalArgumentException("Cannot decrease capacity below the edge's current flow.");
            this.capacity = capacity;
        }

        public int getFlow() {
            return flow;
        }

        public void setFlow(int flow) {
            /* Check that the flow is nonnegative; it's not defined otherwise. */
            if (flow < 0)
                throw new IllegalArgumentException("Flow must be nonnegative.");
            /* Check that the flow doesn't exceed the capacity. */
            if (flow > capacity)
                throw new IllegalArgumentException("Cannot set flow along an edge of capacity " + capacity + " to " + flow);
            this.flow = flow;
        }

        private Edge(T start, T end, int capacity, int flow) {
            this.start = start;
            this.end = end;
            this.capacity = capacity;
            this.flow = flow;
        }

    }

    /* A map from nodes to a secondary map associating other nodes in the graph
     * with the edge connecting the current node to the indicated node.
     */
    private final Map<T, Map<T, Edge<T>>> graph = new HashMap<T, Map<T, Edge<T>>>();

    public boolean addNode(T node) {
        /* If the node already exists, don't do anything. */
        if (graph.containsKey(node))
            return false;
        /* Otherwise, add the node with an empty list of outgoing edges. */
        graph.put(node, new HashMap<T, Edge<T>>());
        return true;
    }

    public Edge<T> addEdge(T start, T dest) {
        /* Confirm both endpoints exist. */
        if (!graph.containsKey(start) || !graph.containsKey(dest))
            throw new NoSuchElementException("Both nodes must be in the graph.");
        /* If this edge doesn't exist, go create it. */
        if (!graph.get(start).containsKey(dest))
            graph.get(start).put(dest, new Edge<T>(start, dest, 0, 0));
        return graph.get(start).get(dest);
    }

    public Edge<T> getEdge(T start, T dest) {
        /* Confirm both endpoints exist. */
        if (!graph.containsKey(start) || !graph.containsKey(dest))
            throw new NoSuchElementException("Both nodes must be in the graph.");
        return graph.get(start).get(dest);
    }

    public void removeEdge(T start, T dest) {
        /* Confirm both endpoints exist. */
        if (!graph.containsKey(start) || !graph.containsKey(dest))
            throw new NoSuchElementException("Both nodes must be in the graph.");
        /* Remove the edge from the graph. */
        graph.get(start).remove(dest);
    }

    public Collection<Edge<T>> edgesFrom(T node) {
        /* Check that the node exists. */
        Map<T, Edge<T>> arcs = graph.get(node);
        if (arcs == null)
            throw new NoSuchElementException("Source node does not exist.");
        return Collections.unmodifiableCollection(arcs.values());
    }

    public Iterator<T> iterator() {
        return graph.keySet().iterator();
    }

    public boolean containsNode(T node) {
        return graph.containsKey(node);
    }

    public int size() {
        return graph.size();
    }
    public boolean isEmpty() {
        return graph.isEmpty();
    }
}

final class ResidualGraph<T> implements Iterable<T> {
    public static final class Edge<T> {
        private final T start;            // The start and end nodes of the edge
        private final T end;
        private final boolean isResidual; // Whether this is a real or residual
                                          // edge
        private Edge<T> reverse;          // The reverse edge for this edge
        private int capacity;             // The capacity of this edge

        public T getStart() {
            return start;
        }

        public T getEnd() {
            return end;
        }

        public int getCapacity() {
            return capacity;
        }

        public void addFlow(int amount) {
            if (amount < 0) {
                reverse.addFlow(-amount);
                return;
            }
            /* Check whether we have enough capacity to support this flow. */
            if (amount > capacity)
                throw new IllegalArgumentException("Cannot push " + amount + " units of flow across edge of capacity " + capacity);
            /* Subtract the appropriate amount from this edge, then add it in
             * to the reverse edge.
             */
            capacity -= amount;
            reverse.capacity += amount;
        }

        public Edge<T> getReverse() {
            return reverse;
        }

        public boolean isOriginal() {
            return !isResidual;
        }

        private void setReverse(Edge<T> reverse) {
            this.reverse = reverse;
        }

        private Edge(T start, T end, int capacity, boolean isOriginal) {
            this.start = start;
            this.end = end;
            this.capacity = capacity;
            this.isResidual = !isOriginal;
        }
    }

    private final Map<T, List<Edge<T>>> graph = new HashMap<T, List<Edge<T>>>();

    public ResidualGraph(FlowNetwork<T> g) {
        for (T node: g)
            graph.put(node, new ArrayList<Edge<T>>());
        for (T node: g) {
            for (FlowNetwork.Edge<T> edge: g.edgesFrom(node)) {
                /* Construct the forward edge, which has capacity equal to the
                 * remaining capacity on this edge.
                 */
                Edge<T> forward = new Edge<T>(edge.getStart(), edge.getEnd(),
                                              edge.getCapacity() - edge.getFlow(),
                                              true);
                /* Construct the reverse edge, whose capacity is the flow on
                 * this edge.
                 */
                Edge<T> reverse = new Edge<T>(edge.getEnd(), edge.getStart(),
                                              edge.getFlow(), false);

                /* Link the two together. */
                forward.setReverse(reverse);
                reverse.setReverse(forward);
                /* Add both to the graph. */
                graph.get(edge.getStart()).add(forward);
                graph.get(edge.getEnd()).add(reverse);
            }
        }
    }

    public Iterator<T> iterator() {
        return graph.keySet().iterator();
    }

    public List<Edge<T>> edgesFrom(T node) {
        /* Look up the list of edges, reporting an error if nothing is found. */
        List<Edge<T>> edges = graph.get(node);
        if (edges == null)
            throw new NoSuchElementException("Node " + node + " does not exist.");
        return Collections.unmodifiableList(edges);
    }
}