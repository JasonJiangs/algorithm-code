package as10;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Stack;

public class TopologicalSort {
//	static String File_Path = "C:/Users/dcdn/Desktop/schedule.txt";
//	static String File_Path = "C:/Users/dcdn/Desktop/four year plan for CPS.txt";
	static String File_Path = "C:/Users/dcdn/Desktop/non-DAG.txt";
	static Graph G;
	static ArrayList<String> courses;

	public static void main(String[] args) throws IOException {
		ArrayList<Integer> resultDFS;
		ArrayList<Integer> resultKahn;
		courses = new ArrayList<String>();
		// obtain all the course and store into courses
		getCourse();
		// create and initialize a graph with size of courses
		G = new Graph(courses.size());
		// to get the connection, we build the graph with prerequisite
		buildGraph();
		resultDFS = G.modifiedDFS(); // store the result of modified DFS algorithm
		resultKahn = G.Kahn(); // store the result of Kahn algorithm
		// check if the result of both algorithm is right
		if(resultDFS.size()!=courses.size()||resultKahn.size()!=courses.size()) {
			System.out.println("IMPOSSIBLE");
		}
		// print the result
		System.out.println("All the courses is listed: ");
		for (int i = 0; i < courses.size(); i++) {
			System.out.print(courses.get(i) + "  ");
		}
		System.out.println("\n");
		// output result of DFS
		System.out.println("The result of topological DFS algorithm is: ");
		for (int i = 0; i < resultDFS.size() - 1; i++) {
			System.out.print(courses.get(resultDFS.get(i)) + " ->");
		}
		System.out.println(courses.get(resultDFS.get(resultDFS.size() - 1)) + "\n");
		// output result of Kahn
		System.out.println("The result of topological Kagn algorithm is: ");
		for (int i = 0; i < resultKahn.size() - 1; i++) {
			System.out.print(courses.get(resultKahn.get(i)) + " ->");
		}
		System.out.println(courses.get(resultKahn.get(resultKahn.size() - 1)) + "\n");
	}
	
	// get the subject that needs to be taken
	public static void getCourse() throws IOException {
		BufferedReader reader;
		try {
			reader = new BufferedReader(new FileReader(File_Path));
			String line = reader.readLine();
			while (line != null) {
				// if the line contain the CLASS which is aim, so
				if (line.contains("CLASS")) {
					String[] words = line.split("\\s+"); // re, split words
					String name = words[words.length - 1];
					courses.add(name); // store course name in courses
				}
				line = reader.readLine(); // repeat in the next line
			}
			reader.close();
		} catch (FileNotFoundException e) {
			System.out.println("An error occured.");
			e.printStackTrace();
		}
	}
	
	public static void buildGraph() throws IOException {
		BufferedReader reader;
		// i is used as index of courses, point the prerequisite of which course
		int i = -1;
		try {
			reader = new BufferedReader(new FileReader(File_Path));
			String line = reader.readLine();
			while (line != null) {
				if (line.contains("CLASS")) {
					i++;
				} else if (line.contains("PREREQ")) {
					// if it's prerequisite, then the course number is added to the graph
					String[] words = line.split("\\s+"); // re, split words
					String name = words[words.length - 1];
					// add i-th course to indexof(name)-th course
					// means course indesof(name) can point to i-th course
					G.addEdge(courses.indexOf(name), i);
				}
				line = reader.readLine();
				// repeat until it ends
			}
			reader.close();
		} catch (FileNotFoundException e) {
			System.out.println("An error occured.");
			e.printStackTrace();
		}
	}
}

class Graph {
	public int V; // number of nodes
	public ArrayList<ArrayList<Integer>> adj; // store classes and their prerequisite
	boolean[] visit; // array that indicates if the node is visited or not
	int[] level; //
	int levelrecord;
	boolean change;

	Graph(int v) { // create a graph and initialization
		V = v;
		visit = new boolean[v];
		level = new int[v];
		adj = new ArrayList<ArrayList<Integer>>(v);
		for (int i = 0; i < v; i++) {
			adj.add(new ArrayList<Integer>()); // adj has v arraylists
			visit[i] = false; // all are not visited yet
			level[i] = 0;
		}
	}

	void addEdge(int v, int u) { // add u-th course to v-th course
		adj.get(v).add(u);
	}

	ArrayList<Integer> modifiedDFS() {// v is index
		levelrecord = 0;
		change = false;
		ArrayList<Integer> result = new ArrayList<Integer>();
		Stack<Integer> stack = new Stack<Integer>();
		for (int i = 0; i < V; i++) {
			if (!visit[i]) {
				dfs(i, visit, stack);
			}
		}
		// pop the stack and collect the result
		while (!stack.empty()) {
			result.add(stack.pop());
		}
		return result;
	}

	void dfs(int v, boolean[] visit, Stack<Integer> stack) {
		// focus on vertex v,
		level[v] = levelrecord;
		visit[v] = true; // change visit to true
		int i;
		java.util.Iterator<Integer> it = adj.get(v).iterator();
		while (it.hasNext()) { // go through all the courses that treat v as prerequisite
			i = it.next();
			if (!visit[i]) { // check if i not visited, recursively do it
				dfs(i, visit, stack);
			} else { // if i visited
				if (level[i] == levelrecord) {// i has next and visited check circle
					return;
				}
			}
		}
		stack.push(v); // after all finish, store the result v, aim to pop in order
		change = false;
		if (change == false) {
			levelrecord++;
			change = true;
		}
	}

	ArrayList<Integer> Kahn() {
		ArrayList<Integer> result = new ArrayList<Integer>();
		int[] indegree = new int[V];
		int resIndex;
		int visitedCounter = 0;
		LinkedList<Integer> queue = new LinkedList<Integer>();
		for (int i = 0; i < V; i++) {
			indegree[i] = 0;
		}
		// calculate every vertex's in-degree
		for (int i = 0; i < V; i++) {
			for (int j = 0; j < adj.get(i).size(); j++) {
				indegree[adj.get(i).get(j)]++;
			}
		}
		// find if any vertex has zero in-degree, put into linked list
		for (int i = 0; i < V; i++) {
			if (indegree[i] == 0) {
				queue.add(i);
			}
		}
		// queue is not empty is it is a DAG
		if (queue.isEmpty() == false) {
			do {
				resIndex = queue.removeFirst();
				visitedCounter++; // count the vertex that has been visited
				result.add(resIndex); // add the index of the vertex to the result
				// go through courses that is free, reduce their in-degree by 1
				for (int i = 0; i < adj.get(resIndex).size(); i++) {
					indegree[adj.get(resIndex).get(i)]--;
					// get courses with in-degree is zero
					if (indegree[adj.get(resIndex).get(i)] == 0) {
						queue.add(adj.get(resIndex).get(i));
					}
				}
			} while (queue.isEmpty() == false); // repeat doing
		}
		if (visitedCounter != V) { // check the correctness
			result.add(-1);
			return result;
		}
		return result;
	}
}