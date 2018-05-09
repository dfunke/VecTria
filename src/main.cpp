//
// Created by funke on 11/23/16.
//

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/point_generators_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

struct tVertexInfo {

    tVertexInfo() : tVertexInfo(-1, -1) {}

    tVertexInfo(const unsigned int &i) : tVertexInfo(i, -1) {}

    tVertexInfo(const unsigned int &i, const unsigned int &p) : id(i), partition(p) {}

    unsigned int id;
    unsigned int partition;
};

struct tSimplexInfo {

    tSimplexInfo() : id(-1) {}

    tSimplexInfo(const unsigned int &i) : id(i) {}

    unsigned int id;
};

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

//2D definitions
using Vb_2 = CGAL::Triangulation_vertex_base_with_info_2<tVertexInfo, K>;
using Cb_2 = CGAL::Triangulation_face_base_with_info_2<tSimplexInfo, K>;
using Tds_2 = CGAL::Triangulation_data_structure_2<Vb_2, Cb_2>;
using DT_2 = CGAL::Delaunay_triangulation_2<K, Tds_2>;
using Point_2 = DT_2::Point;
using Points_2 = std::vector<std::pair<Point_2, tVertexInfo>>;
using Fh_2 = DT_2::Face_handle;
using Vh_2 = DT_2::Vertex_handle;
using Circ_2 = CGAL::Circle_2<K>;
using Box_2 = CGAL::Bbox_2;
using Vec_2 = CGAL::Vector_2<K>;

//3D definitions
using Vb_3 = CGAL::Triangulation_vertex_base_with_info_3<tVertexInfo, K>;
using Cb_3 = CGAL::Triangulation_cell_base_with_info_3<tSimplexInfo, K>;
using Tds_3 = CGAL::Triangulation_data_structure_3<Vb_3, Cb_3>;
using DT_3 = CGAL::Delaunay_triangulation_3<K, Tds_3>;
using Point_3 = DT_3::Point;
using Points_3 = std::vector<std::pair<Point_3, tVertexInfo>>;
using Fh_3 = DT_3::Cell_handle;
using Vh_3 = DT_3::Vertex_handle;
using Circ_3 = CGAL::Sphere_3<K>;
using Box_3 = CGAL::Bbox_3;
using Vec_3 = CGAL::Vector_3<K>;

int main() {

    CGAL::Random_points_in_iso_rectangle_2<Point_2> gen(Point_2(0, 0), Point_2(0.5, 1));
    const unsigned int N = 100;

    //left partition
    Points_2 P1;
    P1.reserve(N);

    unsigned int i = 0;
    for (; i < N; ++i)
        P1.push_back(std::make_pair(*gen++, tVertexInfo(i, 0)));
    DT_2 T1(P1.begin(), P1.end());

    std::cout << "Left partition: " << T1.number_of_vertices() << " points, " << T1.number_of_faces() << " faces"
              << std::endl;

    //right partition
    Points_2 P2;
    P2.reserve(N);

    for (; i < 2 * N; ++i)
        P2.push_back(std::make_pair((*gen++) + Vec_2(0.5, 0), tVertexInfo(i, 1)));
    DT_2 T2(P2.begin(), P2.end());

    std::cout << "Left partition: " << T2.number_of_vertices() << " points, " << T2.number_of_faces() << " faces"
              << std::endl;

    return EXIT_SUCCESS;

}