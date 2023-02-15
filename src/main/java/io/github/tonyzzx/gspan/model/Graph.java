package io.github.tonyzzx.gspan.model;

import java.io.BufferedReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class Graph extends ArrayList<Vertex> {
    private static final long serialVersionUID = 1L;
	private static Set<Integer> boringEdges = new HashSet<>(); // new HashSet<>(Arrays.asList(1, 4, 7)); // "def" and "recv", "rep"
	
	
	public static long mustIncludeNode = -1; // init in Main.main(String[])
	public static List<Integer> boringNodes; // init in Main.main(String[])
	
    public int edge_size = 0;
    boolean directed = false;
    
    
    public int ID;
    
    public char label;
    public int quantity = 0; // In this formulation of the problem, each graph represents many clones.

    public Graph() {
    }

    public Graph(boolean directed) {
        this.directed = directed;
    }

    void buildEdge() {
        String buf;
        NavigableMap<String, Integer> tmp = new TreeMap<>();

        int id = 0;
        for (int from = 0; from < size(); ++from) {
            for (Edge it : this.get(from).edge) {
                if (directed || from <= it.to)
                    buf = from + " " + it.to + " " + it.eLabel;
                else
                    buf = it.to + " " + from + " " + it.eLabel;

                // Assign unique id's for the edges.
                if (tmp.get(buf) == null) {
                    it.id = id;
                    tmp.put(buf, id);
                    ++id;
                } else {
                    it.id = tmp.get(buf);
                }
            }
        }

        edge_size = id;
    }

    public BufferedReader read(BufferedReader is) throws IOException {
        ArrayList<String> result = new ArrayList<>();
        String line;

        clear();

        while ((line = is.readLine()) != null) {
            result.clear();
            String[] splitRead = line.split(" ");
            Collections.addAll(result, splitRead);

            if (!result.isEmpty()) {
            	if (result.get(0).equals("-")) { // delimiter
//            		if (!this.isEmpty()) {
//                        break;
//                    } else {
//                    	throw new RuntimeException("Unexpected file format");
//                    }
            		break;
            	} else if (result.get(0).equals("t")) {
            		
            		String supposedID = result.get(2);
            		this.ID = Integer.parseInt(supposedID);
                    
                    char label = result.get(3).charAt(0);
                    // 'U' for unlabeled
                    // 'M' for misuse
                    // 'C' for correct
                    assert label == 'U' || label == '+' || label == '-';
                    this.label = label;
                    
                    int howMany = Integer.parseInt(result.get(4));
                    this.quantity = howMany;
                    assert this.quantity >= 1;
                    
                } else if (result.get(0).equals("v") && result.size() >= 3) {
                    // int id = Integer.parseInt(result.get(1));
                    Vertex vex = new Vertex();
                    vex.label = Integer.parseInt(result.get(2));
                    this.add(vex);
                    
                } else if (result.get(0).equals("e") && result.size() >= 4) {
                    int from = Integer.parseInt(result.get(1));
                    int to = Integer.parseInt(result.get(2));
                    int eLabel = Integer.parseInt(result.get(3));

                    if (this.size() <= from || this.size() <= to) {
                        System.out.println("Format Error:  define vertex lists before edges");
                        return null;
                    }

                    this.get(from).push(from, to, eLabel);

                    if (!directed) {
                        this.get(to).push(to, from, eLabel);
                    }
                }
            }
        }

        buildEdge();

        return is;
    }

    public void write(FileWriter os) throws IOException {
        String buf;
        // Sort the result of edges.
        NavigableSet<String> tmp = new TreeSet<>((o1, o2) -> {
            String[] split1 = o1.split(" ");
            String[] split2 = o2.split(" ");
            if (Integer.parseInt(split1[0]) == Integer.parseInt(split2[0])) {
                if (Integer.parseInt(split1[1]) == Integer.parseInt(split2[1])) {
                    return Integer.parseInt(split1[2]) - Integer.parseInt(split2[2]);
                } else {
                    return Integer.parseInt(split1[1]) - Integer.parseInt(split2[1]);
                }
            } else {
                return Integer.parseInt(split1[0]) - Integer.parseInt(split2[0]);
            }
        });

        for (int from = 0; from < size(); ++from) {
            os.write("v " + from + " " + this.get(from).label + System.getProperty("line.separator"));

            for (Edge it : this.get(from).edge) {
                if (directed || from <= it.to) {
                    buf = from + " " + it.to + " " + it.eLabel;
                } else {
                    buf = it.to + " " + from + " " + it.eLabel;
                }
                tmp.add(buf);
            }
        }

        for (String it : tmp) {
            os.write("e " + it + System.getProperty("line.separator"));
        }

        os.flush();
    }
    
    public boolean allEdgesAreBoring() {
    	if (size() <= 1) { // single node
    		return false;
    	}
    	
    	boolean result = true;
    	for (int from = 0; from < size(); ++from) {

            for (Edge it : this.get(from).edge) {
          
                if (!Graph.boringEdges.contains(it.eLabel)) {
                	result = false;
                }
            
            }
        }
    	return result;
    }
    
    public boolean allNodesAreBoring() {
    	if (size() == 1) {
    		return Graph.boringNodes.contains(this.get(0).label);
    	}

    	
    	boolean result = true;
    	for (int from = 0; from < size(); ++from) {

            if (!Graph.boringNodes.contains(this.get(from).label)) {
            	result = false;
            }
                    
        }
    	return result;
    }
    
    public boolean containsMustIncludeNode() {
    	if (mustIncludeNode == -1) {
    		return true;
    	}

    	for (int from = 0; from < size(); ++from) {
//    		System.out.println("matching... mustIncludeNode = " + mustIncludeNode + " but label is  " + this.get(from).label );
    		if (this.get(from).label == (int)mustIncludeNode) {
//    			System.out.println("matched mustIncludeNode");
    			return true;
    		}
    	}
    	
    	return false;
    }

    public void check() {
        /*
         * Check all indices
         */
        for (int from = 0; from < size(); ++from) {
            System.out.println(
                    "check vertex " + from + ", label " + this.get(from).label + System.getProperty("line.separator"));

            for (Edge it : this.get(from).edge) {
                System.out.println("   check edge from " + it.from + " to " + it.to + ", label " + it.eLabel
                        + System.getProperty("line.separator"));
                assert (it.from >= 0 && it.from < size());
                assert (it.to >= 0 && it.to < size());
            }
        }
    }

    void resize(int size) {
        while (this.size() < size) {
            this.add(new Vertex());
        }
    }
}
