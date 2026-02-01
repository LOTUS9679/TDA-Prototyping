#include <iostream>
#include <vector>
#include <algorithm>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

typedef CGAL::Simple_cartesian<double> K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef K::Point_3 Point;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::edge_descriptor edge_descriptor;

// Helper: Find the sink (Lake Bottom) for any vertex
vertex_descriptor find_sink(vertex_descriptor start, Mesh& mesh) {
    vertex_descriptor current = start;
    while(true) {
        vertex_descriptor lowest = current;
        double min_h = mesh.point(current).z();
        for(auto neighbor : CGAL::vertices_around_target(mesh.halfedge(current), mesh)) {
            if(mesh.point(neighbor).z() < min_h) {
                min_h = mesh.point(neighbor).z();
                lowest = neighbor;
            }
        }
        if(lowest == current) return current;
        current = lowest;
    }
}

int main()
{
    Mesh mesh;
    int size = 5; 
    std::vector<std::vector<vertex_descriptor> > verts(size, std::vector<vertex_descriptor>(size));

    // 1. Build Grid
    for(int i=0; i<size; ++i) 
        for(int j=0; j<size; ++j) 
            verts[i][j] = mesh.add_vertex(Point(i, j, 0));
            
    for(int i=0; i<size-1; ++i) 
        for(int j=0; j<size-1; ++j) 
            mesh.add_face(verts[i][j], verts[i+1][j], verts[i+1][j+1], verts[i][j+1]);

    // 2. Hardcode Heights (Same as before)
    std::vector<int> heights = {
        3, 1, 4, 0, 8,
        4, 4, 6, 8, 4,
        5, 5, 0, 7, 6,
        9, 4, 5, 4, 4,
        3, 9, 0, 8, 9
    };
    int idx = 0;
    for(vertex_descriptor v : vertices(mesh)) mesh.point(v) = Point(0, 0, heights[idx++]);

    std::cout << "--- SEARCHING FOR SADDLES (Border Crossings) ---" << std::endl;

    // 3. Scan all Edges to find "Borders"
    for(edge_descriptor e : edges(mesh)) {
        vertex_descriptor u = source(e, mesh);
        vertex_descriptor v = target(e, mesh);

        // Find which Lake each endpoint belongs to
        vertex_descriptor sink_u = find_sink(u, mesh);
        vertex_descriptor sink_v = find_sink(v, mesh);

        // If they belong to DIFFERENT Lakes, this edge is a SADDLE
        if (sink_u != sink_v) {
            double height_u = mesh.point(u).z();
            double height_v = mesh.point(v).z();
            
            // The "Height" of the pass is the higher of the two points
            double saddle_height = std::max(height_u, height_v);
            
            // Calculate Persistence: Saddle Height - (Higher Sink Height)
            // (Note: In standard persistence, the "younger" feature dies.
            //  The younger feature is usually the one with the HIGHER birth/sink height).
            double sink_h_u = mesh.point(sink_u).z();
            double sink_h_v = mesh.point(sink_v).z();
            double higher_sink_h = std::max(sink_h_u, sink_h_v);
            
            double persistence = saddle_height - higher_sink_h;

            std::cout << "  SADDLE FOUND between Lake " << sink_u << " & Lake " << sink_v << std::endl;
            std::cout << "    Saddle Height: " << saddle_height << std::endl;
            std::cout << "    Persistence: " << persistence << std::endl;
            
            if (persistence < 3.0) {
                 std::cout << " VERDICT: NOISE (Merge lakes)" << std::endl;
            } else {
                 std::cout << " VERDICT: FEATURE (Keep separate)" << std::endl;
            }
            std::cout << "-----------------------------------" << std::endl;
        }
    }
    
    return 0;
}
