//
// Created by funke on 11/23/16.
//

#include <random>
#include <functional>
#include <chrono>

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

#ifdef HAS_VTUNE
#include <ittnotify.h>
#endif

#include "GeometryStructures.h"
#include "Predicates.h"

#define SEED 1986

template<tDimType D, typename Precision, class PointArray>
void generatePoints(PointArray &pa, const tIndexType &n) {

    std::mt19937 generator(SEED);
    std::uniform_real_distribution<Precision> distribution;
    auto rand = std::bind(distribution, generator);

    pa.ensure(n);
    for (tIndexType i = 0; i < n; ++i) {
        auto p = pa.get(i);

        for (tDimType d = 0; d < D; ++d) {
            p[d] = rand();
        }
    }
}

template<tDimType D>
struct Triangulator;

template<>
struct Triangulator<2> {

    template<class SimplexArray, class PointArray>
    SimplexArray triangulate(const PointArray &pa) {

        using K = CGAL::Exact_predicates_inexact_constructions_kernel;

//2D definitions
        using Vb_2 = CGAL::Triangulation_vertex_base_with_info_2<tIndexType, K>;
        using Cb_2 = CGAL::Triangulation_face_base_with_info_2<tIndexType, K>;
        using Tds_2 = CGAL::Triangulation_data_structure_2<Vb_2, Cb_2>;
        using DT_2 = CGAL::Delaunay_triangulation_2<K, Tds_2>;
        using Point_2 = DT_2::Point;
        using Points_2 = std::vector<std::pair<Point_2, tIndexType>>;
//        using Fh_2 = DT_2::Face_handle;
//        using Vh_2 = DT_2::Vertex_handle;
//        using Circ_2 = CGAL::Circle_2<K>;
//        using Box_2 = CGAL::Bbox_2;
//        using Vec_2 = CGAL::Vector_2<K>;

        Points_2 P;
        P.reserve(pa.size());

        for (unsigned int i = 0; i < pa.size(); ++i) {
            auto p = pa.get(i);
            Point_2 p_cgal(p[0], p[1]);
            P.push_back(std::make_pair(p_cgal, i));
        }

        DT_2 T(P.begin(), P.end());

        tIndexType simplexId = 0;
        for (auto it = T.finite_faces_begin(); it != T.finite_faces_end(); ++it) {
            it->info() = simplexId++;
        }

        SimplexArray simplices;
        simplices.ensure(simplexId - 1);

        for (auto it = T.finite_faces_begin(); it != T.finite_faces_end(); ++it) {
            tIndexType id = it->info();
            auto s = simplices.get(id);

            s.vertex(0) = it->vertex(0)->info();
            s.vertex(1) = it->vertex(1)->info();
            s.vertex(2) = it->vertex(2)->info();

            s.neighbor(0) = T.is_infinite(it->neighbor(0)) ? INF : it->neighbor(0)->info();
            s.neighbor(1) = T.is_infinite(it->neighbor(1)) ? INF : it->neighbor(1)->info();
            s.neighbor(2) = T.is_infinite(it->neighbor(2)) ? INF : it->neighbor(2)->info();
        }

        return simplices;

    }
};

template<>
struct Triangulator<3> {

    template<class SimplexArray, class PointArray>
    SimplexArray triangulate(const PointArray &pa) {

        using K = CGAL::Exact_predicates_inexact_constructions_kernel;

//3D definitions
        using Vb_3 = CGAL::Triangulation_vertex_base_with_info_3<tIndexType, K>;
        using Cb_3 = CGAL::Triangulation_cell_base_with_info_3<tIndexType, K>;
        using Tds_3 = CGAL::Triangulation_data_structure_3<Vb_3, Cb_3>;
        using DT_3 = CGAL::Delaunay_triangulation_3<K, Tds_3>;
        using Point_3 = DT_3::Point;
        using Points_3 = std::vector<std::pair<Point_3, tIndexType>>;
//        using Fh_3 = DT_3::Cell_handle;
//        using Vh_3 = DT_3::Vertex_handle;
//        using Circ_3 = CGAL::Sphere_3<K>;
//        using Box_3 = CGAL::Bbox_3;
//        using Vec_3 = CGAL::Vector_3<K>;

        Points_3 P;
        P.reserve(pa.size());

        for (unsigned int i = 0; i < pa.size(); ++i) {
            auto p = pa.get(i);
            Point_3 p_cgal(p[0], p[1], p[2]);
            P.push_back(std::make_pair(p_cgal, i));
        }

        DT_3 T(P.begin(), P.end());

        tIndexType simplexId = 0;
        for (auto it = T.finite_cells_begin(); it != T.finite_cells_end(); ++it) {
            it->info() = simplexId++;
        }

        SimplexArray simplices;
        simplices.ensure(simplexId - 1);

        for (auto it = T.finite_cells_begin(); it != T.finite_cells_end(); ++it) {
            tIndexType id = it->info();
            auto s = simplices.get(id);

            s.vertex(0) = it->vertex(0)->info();
            s.vertex(1) = it->vertex(1)->info();
            s.vertex(2) = it->vertex(2)->info();
            s.vertex(3) = it->vertex(3)->info();

            s.neighbor(0) = T.is_infinite(it->neighbor(0)) ? INF : it->neighbor(0)->info();
            s.neighbor(1) = T.is_infinite(it->neighbor(1)) ? INF : it->neighbor(1)->info();
            s.neighbor(2) = T.is_infinite(it->neighbor(2)) ? INF : it->neighbor(2)->info();
            s.neighbor(3) = T.is_infinite(it->neighbor(3)) ? INF : it->neighbor(3)->info();
        }

        return simplices;

    }
};

template<tDimType D, typename Precision>
struct Checker;

template<typename Precision>
struct Checker<3, Precision> {

    static constexpr tDimType D = 3;

    template<class SimplexArray, class PointArray>
    bool check(const SimplexArray &simplices, const PointArray &points) {

        bool valid = true;

        for (tIndexType i = 0; i < simplices.size(); ++i) {
            auto s = simplices.get(i);

            for (tDimType n = 0; n < D + 1; ++n) {
                if (s.neighbor(n) != INF) {
                    auto sn = simplices.get(n);

                    for (tDimType j = 0; j < D + 1; ++j) {
                        Precision det = insphere_fast<D, Precision>(s.vertex(0), s.vertex(1), s.vertex(2), s.vertex(3),
                                                                    sn.vertex(j), points);

                        if (det < 0) {
                            valid = false;
                        }
                    }
                }
            }
        }

        return valid;
    }
};

#define D 3
#define Precision double
#define N 1e6

int main() {

#ifdef HAS_VTUNE
    __itt_pause();
#endif

    PointArray<D, Precision, MemoryLayoutAoA> points_aoa;
    generatePoints<D, Precision>(points_aoa, N);

    PointArray<D, Precision, MemoryLayoutPA> points_pa;
    generatePoints<D, Precision>(points_pa, N);

    Triangulator<D> triangulator;

    auto simplices_aoa = triangulator.triangulate<SimplexArray<D, MemoryLayoutAoA>>(points_aoa);
    auto simplices_pa = triangulator.triangulate<SimplexArray<D, MemoryLayoutPA>>(points_pa);

    bool valid;
    Checker<D, Precision> checker;

#ifdef HAS_VTUNE
    __itt_resume();
#endif

    auto t1 = std::chrono::high_resolution_clock::now();
    valid = checker.check(simplices_aoa, points_aoa);
    auto t2 = std::chrono::high_resolution_clock::now();

    std::cout << "sAOA/pAOA valid: " << valid << " time "
              << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    valid = checker.check(simplices_aoa, points_pa);
    t2 = std::chrono::high_resolution_clock::now();

    std::cout << "sAOA/pPA valid: " << valid << " time "
              << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    valid = checker.check(simplices_pa, points_aoa);
    t2 = std::chrono::high_resolution_clock::now();

    std::cout << "sPA/pAOA valid: " << valid << " time "
              << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    valid = checker.check(simplices_pa, points_pa);
    t2 = std::chrono::high_resolution_clock::now();

    std::cout << "sPA/pPA valid: " << valid << " time "
              << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() << std::endl;

    return 0;

}
