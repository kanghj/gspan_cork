package io.github.tonyzzx.gspan.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.Vector;

public class EfficientHistory {

	private Set<Integer> edge;
    private Set<Integer> vertex;
    
    public List<Edge> ordering;

    public EfficientHistory(Graph g, PDFS p) {
//    	System.out.println("\t\tConstructing EfficientHistory : " + g.edge_size + ", " + g.size());
        edge = new HashSet<>(g.edge_size + 1, 1.0f);
        vertex = new HashSet<>(g.size() + 1, 1.0f);
        ordering = new ArrayList<>(g.edge_size + 1);
        
        build(g, p);
    }

    private void build(Graph graph, PDFS e) {
        // first build history
    	ordering.clear();
        edge.clear();
//        edge.setSize(graph.edge_size);
        vertex.clear();
//        vertex.setSize(graph.size());


        if (e != null) {
        	ordering.add(e.edge);
            edge.add(e.edge.id);
  
            
            vertex.add(e.edge.from);
            vertex.add(e.edge.to);

            for (PDFS p = e.prev; p != null; p = p.prev) {
            	ordering.add(p.edge); // this line eats 8% of overall instructions(!)
                edge.add(p.edge.id);
      
                
                
                vertex.add(p.edge.from);
                vertex.add(p.edge.to);
            }
            Collections.reverse(ordering);
        }
    }

    public boolean hasEdge(int id) {
        return edge.contains(id);
    }

    public boolean hasVertex(int id) {
        return vertex.contains(id);
    }
}
