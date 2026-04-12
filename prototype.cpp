// CGAL Tools
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/IO/polygon_mesh_io.h>

// Gudhi Tools
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>

// Standard C++ Tools
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <CGAL/Random.h>

// The kernel and mesh types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Mesh;
typedef K::Point_3 Point;
typedef Mesh::Vertex_index Vertex_index;
typedef Mesh::Edge_index Edge_index;

// For CGAL's benchmark method
typedef boost::property_map<Mesh, CGAL::edge_is_feature_t>::type EIFMap;

// ============================================================
// BENCHMARK FUNCTION
// Tests CGAL built-in method and counts sharp edges
// ============================================================
int run_cgal_benchmark(Mesh& mesh) {
    EIFMap eif = get(CGAL::edge_is_feature, mesh);
    
    // 45 = threshold angle in degrees
    CGAL::Polygon_mesh_processing::detect_sharp_edges(mesh, 45, eif);
    
    int sharp_count = 0;
    for(auto e : mesh.edges()) {
        if(get(eif, e)) {
            sharp_count++;
        }
    }
    return sharp_count;
}

// ============================================================
// CALCULATE DIHEDRAL ANGLE FILTRATION
// Dihedral angle = angle between two faces sharing an edge
// Sharp edge = large dihedral angle deviation from flat
// ============================================================
double calculate_filtration(Mesh& mesh, Edge_index e) {
    auto h1 = mesh.halfedge(e, 0);
    auto h2 = mesh.halfedge(e, 1);

    // Boundary edges are not sharp features in this context
    if(mesh.is_border(h1) || mesh.is_border(h2)) {
        return 0.0;
    }

    auto face1 = mesh.face(h1);
    auto face2 = mesh.face(h2);

    if(face1 == Mesh::null_face() || face2 == Mesh::null_face()) {
        return 0.0;
    }

    // Get vertices of face 1
    auto f1_h = mesh.halfedge(face1);
    auto v1_0 = mesh.point(mesh.source(f1_h));
    auto v1_1 = mesh.point(mesh.target(f1_h));
    auto v1_2 = mesh.point(mesh.target(mesh.next(f1_h)));

    // Get vertices of face 2
    auto f2_h = mesh.halfedge(face2);
    auto v2_0 = mesh.point(mesh.source(f2_h));
    auto v2_1 = mesh.point(mesh.target(f2_h));
    auto v2_2 = mesh.point(mesh.target(mesh.next(f2_h)));

    // Calculate normals
    auto edge1_f1 = v1_1 - v1_0;
    auto edge2_f1 = v1_2 - v1_0;
    auto normal1 = CGAL::cross_product(edge1_f1, edge2_f1);

    auto edge1_f2 = v2_1 - v2_0;
    auto edge2_f2 = v2_2 - v2_0;
    auto normal2 = CGAL::cross_product(edge1_f2, edge2_f2);

    double len1 = std::sqrt(CGAL::to_double(normal1.squared_length()));
    double len2 = std::sqrt(CGAL::to_double(normal2.squared_length()));

    if(len1 == 0.0 || len2 == 0.0) return 0.0;

    double dot = CGAL::to_double(normal1 * normal2);
    double cos_angle = dot / (len1 * len2);
    cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
    
    double angle = std::acos(cos_angle);

    // Return deviation from flat (M_PI = flat surface)
    return std::abs(M_PI - angle);
}
// ============================================================
// RUN TDA PIPELINE
// ============================================================
int run_tda_pipeline(Mesh& mesh, double angle_threshold) {
    std::cout << "  Building Simplex Tree..." << std::endl;
    Gudhi::Simplex_tree<> st;

    for(auto v : mesh.vertices()) {
        st.insert_simplex({(int)v}, 0.0);
    }

    int preserved_edges = 0;

    for(auto e : mesh.edges()) {
        Vertex_index v0 = mesh.vertex(e, 0);
        Vertex_index v1 = mesh.vertex(e, 1);
        
        // This is the value Gudhi uses to build the topology
        double filtration = calculate_filtration(mesh, e);
        st.insert_simplex({(int)v0, (int)v1}, filtration);

        // Calculate the raw angle back from the filtration value
        double actual_angle = M_PI - filtration;
        
        // Count it ONLY if the angle is sharper than our threshold
        if(actual_angle >= angle_threshold) {
            preserved_edges++;
        }
    }

    st.make_filtration_non_decreasing();

    std::cout << "  Computing persistence..." << std::endl;
    typedef Gudhi::persistent_cohomology::Field_Zp Field;
    typedef Gudhi::persistent_cohomology::Persistent_cohomology<Gudhi::Simplex_tree<>, Field> Persistence;

    Persistence pcoh(st);
    pcoh.init_coefficients(2);
    pcoh.compute_persistent_cohomology();

    int tda_features = 0;
    auto intervals = pcoh.get_persistent_pairs();

    // Low threshold to successfully group the connected topological features
    double persistence_threshold = 0.05; 

    for(auto& interval : intervals) {
        auto birth_simplex = get<0>(interval);
        auto death_simplex = get<1>(interval);

        if(death_simplex == st.null_simplex()) continue;

        double birth = st.filtration(birth_simplex);
        double death = st.filtration(death_simplex);
        double lifespan = std::abs(death - birth);

        if(lifespan >= persistence_threshold) {
            tda_features++;
        }
    }

    std::cout << "  -> TDA grouped these edges into " << tda_features << " macroscopic topological features." << std::endl;

    return preserved_edges;
}
// ============================================================
// MAIN FUNCTION
// ============================================================
int main() {
    
        Mesh mesh;

    // STAGE 1: LOAD SEBASTIEN'S FILE
    std::cout << "Stage 1: Loading mentor's angle_test.off file..." << std::endl;
    if(!CGAL::IO::read_polygon_mesh("angle_test.off", mesh)) {
        std::cerr << "CRITICAL ERROR: Could not read angle_test.off!" << std::endl;
        return 1;
    }
    std::cout << "Mesh loaded successfully!" << std::endl;

    // STAGE 2: RUN THE CGAL BENCHMARK
    std::cout << "Stage 2: Running standard CGAL benchmark..." << std::endl;
    int cgal_sharp_count = run_cgal_benchmark(mesh);
    std::cout << "CGAL detect_sharp_edges found : " << cgal_sharp_count << " sharp edges" << std::endl;
    std::cout << "====================================" << std::endl;

    // STAGE 3: RUN TDA PIPELINE
    std::cout << "Stage 3: Running Geometry-Aware TDA pipeline..." << std::endl;
    
    // We set the threshold to 1.2 (approx 68 degrees) to safely filter out the shallow curved noise
    double threshold = 1.2; 
    int tda_sharp_count = run_tda_pipeline(mesh, threshold);
    std::cout << "TDA pipeline complete!" << std::endl;

    
  /*
   Mesh mesh;

    // STAGE 1: BUILD THE 20x20 FLAT GRID
    std::cout << "Stage 1: Building 20x20 flat grid..." << std::endl;
    int size = 20;
    std::vector<std::vector<Vertex_index> > verts(size, std::vector<Vertex_index>(size));

    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {
            verts[i][j] = mesh.add_vertex(Point(i, j, 0));
        }
    }

    std::cout << "Stage 2: Connecting faces..." << std::endl;
    for(int i = 0; i < size - 1; i++) {
        for(int j = 0; j < size - 1; j++) {
            mesh.add_face(verts[i][j], verts[i+1][j], verts[i+1][j+1], verts[i][j+1]);
        }
    }

    // STAGE 3: ADD RANDOM NOISE (-0.4 to +0.4)
    std::cout << "Stage 3: Adding random Z-noise..." << std::endl;
    CGAL::Random random(42);
    for(auto v : mesh.vertices()) {
        Point p = mesh.point(v);
        mesh.point(v) = Point(p.x(), p.y(), random.get_double(-0.4, 0.4));
    }
    

    // STAGE 4: RUN THE CGAL BENCHMARK
    std::cout << "Stage 4: Running standard CGAL benchmark..." << std::endl;
    int cgal_sharp_count = run_cgal_benchmark(mesh);
    std::cout << "CGAL detect_sharp_edges found : " << cgal_sharp_count << " false sharp edges" << std::endl;
    std::cout << "====================================" << std::endl;

    // STAGE 5: RUN TDA PIPELINE
    std::cout << "Stage 5: Running Geometry-Aware TDA pipeline..." << std::endl;
    
    // We keep the threshold at 1.2 (approx 68 degrees) to filter the shallow noise bumps
    double threshold = 1.2; 
    int tda_sharp_count = run_tda_pipeline(mesh, threshold);
    std::cout << "TDA pipeline complete!" << std::endl;

*/

    // STAGE 4: FINAL COMPARISON OUTPUT
    std::cout << std::endl;
    std::cout << "========================================" << std::endl;
    std::cout << "       FINAL COMPARISON RESULTS         " << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "CGAL detect_sharp_edges (old method) : " << cgal_sharp_count << " edges detected" << std::endl;
    std::cout << "TDA Persistence method  (new method) : " << tda_sharp_count << " real edges preserved" << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    int noise_removed = cgal_sharp_count - tda_sharp_count;
    std::cout << "False positives filtered by TDA      : " << noise_removed << " edges removed!" << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "Total Vertices : " << mesh.number_of_vertices() << std::endl;
    std::cout << "Total Edges    : " << mesh.number_of_edges() << std::endl;
    std::cout << "Total Faces    : " << mesh.number_of_faces() << std::endl;
    std::cout << "Benchmark run successfully!" << std::endl; 
    
    return 0;
}