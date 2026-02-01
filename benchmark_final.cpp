#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Random.h>
#include <boost/property_map/property_map.hpp> // Added for safety
#include <iostream>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef K::Point_3 Point;

int main() {
    Mesh mesh;
    int size = 20; 
    std::vector<std::vector<Mesh::Vertex_index>> verts(size, std::vector<Mesh::Vertex_index>(size));

    // 1. Create a Flat Grid
    for(int i=0; i<size; ++i) {
        for(int j=0; j<size; ++j) {
            verts[i][j] = mesh.add_vertex(Point(i, j, 0));
        }
    }
    for(int i=0; i<size-1; ++i) {
        for(int j=0; j<size-1; ++j) {
            mesh.add_face(verts[i][j], verts[i+1][j], verts[i+1][j+1], verts[i][j+1]);
        }
    }

    // 2. Add Random Noise (Amplitude 0.4)
    CGAL::Random random(42);
    for(auto v : mesh.vertices()) {
        Point p = mesh.point(v);
        mesh.point(v) = Point(p.x(), p.y(), random.get_double(-0.4, 0.4));
    }

    // 3. Detect Features
    typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;
    EIFMap eif = get(CGAL::edge_is_feature, mesh);
    int sharp_count = 0;

    // Run detection with 45 degree threshold
    CGAL::Polygon_mesh_processing::detect_sharp_edges(mesh, 45, eif);

    for(auto e : mesh.edges()) {
        if(get(eif, e)) sharp_count++;
    }

    std::cout << "Detected " << sharp_count << " sharp edges on noisy plane." << std::endl;
    return 0;
}
