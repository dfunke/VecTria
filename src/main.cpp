//
// Created by funke on 11/23/16.
//

#include <vector>
#include <unordered_set>
#include <algorithm>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/Dimension.h>

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

using sVh_2 = std::unordered_set<Vh_2>;
using sFh_2 = std::unordered_set<Fh_2>;

unsigned int pointId = 0;
unsigned int partitionId = 0;

/* Tests whehter a sphere intersects with the box */
template <typename BOX, typename SPHERE>
static bool boxIntersects(const BOX &box, const SPHERE &sphere) {
    const auto D = CGAL::Ambient_dimension<BOX>::value;
    auto r2 = sphere.squared_radius();
    auto r = std::sqrt(r2);

    double dist = 0;
    for (unsigned short d = 0; d < D; ++d) {
        double e = std::max(box.min(d) - sphere.center()[d], (double)0) +
                   std::max(sphere.center()[d] - box.max(d), (double)0);

        if (r < e) return false;
        dist += e * e;
    }

    if (dist <= r2) return true;

    return false;
}

/* tests whether sphere is FULLY contained in box */
template <typename BOX, typename SPHERE>
static bool boxContains(const BOX &box, const SPHERE &sphere) {
    const auto D = CGAL::Ambient_dimension<BOX>::value;
    auto r2 = sphere.squared_radius();

    // check whether center is in box
    for (unsigned short d = 0; d < D; ++d) {
        if (!(box.min(d) <= sphere.center()[d] &&
              sphere.center()[d] <= box.max(d))) {
            return false;
        }
    }

    // the center of the sphere is within the box
    for (unsigned short d = 0; d < D; ++d) {
        // project p to the boundary of box in dimension d closest to center of
        // the sphere
        double dist =
                sphere.center()[d] - (sphere.center()[d] < (box.max(d) + box.min(d)) / 2
                                      ? box.min(d)
                                      : box.max(d));
        dist = dist * dist;

        if (dist < r2) return false;
    }

    return true;
}

template<class Gen>
auto generate(Gen & gen, const unsigned int N, const Vec_2 & offset){

    Points_2 P;
    P.reserve(N);

    for (unsigned int i = 0; i < N; ++i)
        P.push_back(std::make_pair(*gen++ + offset, tVertexInfo(pointId++, partitionId)));
    DT_2 T(P.begin(), P.end());

    std::cout << "partition = " << partitionId
              << ", points = " << T.number_of_vertices()
              << ", faces = " << T.number_of_faces()
              << std::endl;

    partitionId++;

    return std::make_tuple(P, T);
}

auto bounds(const Points_2 &P){
    auto x = std::minmax_element(P.begin(), P.end(), [](const auto & a, const auto & b) -> bool { return a.first.x() < b.first.x();});
    auto y = std::minmax_element(P.begin(), P.end(), [](const auto & a, const auto & b) -> bool { return a.first.y() < b.first.y();});

    return Box_2(x.first->first.x(), y.first->first.y(), x.second->first.x(), y.second->first.y());
}

auto border(DT_2 T, const Box_2 & b){
    sVh_2 Vb;
    sFh_2 Fb;

    for(auto f = T.finite_faces_begin(); f != T.finite_faces_end(); ++f){
        Circ_2 c(f->vertex(0)->point(),
                 f->vertex(1)->point(),
                 f->vertex(2)->point());

        if(boxIntersects(b, c)){
            Fb.insert(f);

            Vb.insert(f->vertex(0));
            Vb.insert(f->vertex(1));
            Vb.insert(f->vertex(2));
        }
    }

    std::cout << "partition = " << T.finite_vertices_begin()->info().partition
              << ", border points = " << Vb.size()
              << ", border faces = " << Fb.size()
              << std::endl;

    return std::make_tuple(Vb, Fb);
}

auto triangulateBorder(const sVh_2 & P1, const sVh_2 & P2){


    auto transform = [] (const Vh_2 & v) {
        return std::make_pair(v->point(), v->info());
    };

    DT_2 T;
    T.insert(boost::make_transform_iterator(P1.begin(), transform), boost::make_transform_iterator(P1.end(), transform));
    T.insert(boost::make_transform_iterator(P2.begin(), transform), boost::make_transform_iterator(P2.end(), transform));

    std::cout << "border"
              << ", points = " << T.number_of_vertices()
              << ", faces = " << T.number_of_faces()
              << std::endl;

    return T;
}

auto mergeTriangulations(const DT_2 & T1, const sFh_2 & Fb1,
                         const DT_2 & T2, const sFh_2 & Fb2,
                         const DT_2 & Tb){

    DT_2 T;

    for(auto it = T1.finite_faces_begin(); it != T1.finite_faces_end(); ++it){
        if(!Fb1.count(it))
            T.create_face(it);
    }
    for(auto it = T2.finite_faces_begin(); it != T2.finite_faces_end(); ++it){
        if(!Fb2.count(it))
            T.create_face(it);
    }

    for(auto it = Tb.finite_faces_begin(); it != Tb.finite_faces_end(); ++it){
        if(it->vertex(0)->info().partition != it->vertex(0)->info().partition
                || it->vertex(1)->info().partition != it->vertex(2)->info().partition
                || it->vertex(0)->info().partition != it->vertex(2)->info().partition){

            T.create_face(it);
        }
    }

    std::cout << "merged"
              << ", points = " << T.number_of_vertices()
              << ", faces = " << T.number_of_faces()
              << std::endl;

    return T;
}

int main() {

    CGAL::Random_points_in_iso_rectangle_2<Point_2> gen(Point_2(0, 0), Point_2(0.5, 1));
    const unsigned int N = 100;

    //left partition
    auto [P1, T1] = generate(gen, N, Vec_2(0,0));

    //right partition
    auto [P2, T2] = generate(gen, N, Vec_2(0.5,0));

    auto [Pb1, Sb1] = border(T1, bounds(P2));
    auto [Pb2, Sb2] = border(T2, bounds(P1));

    auto Tb = triangulateBorder(Pb1, Pb2);
    auto T = mergeTriangulations(T1, Sb1, T2, Sb2, Tb);

    return EXIT_SUCCESS;

}