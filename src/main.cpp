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

#ifdef HAS_ADVISOR
#include <advisor-annotate.h>
#endif

#include "GeometryStructures.h"
#include "Predicates.h"

#define SEED 1986

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

//2D definitions
using Vb_2 = CGAL::Triangulation_vertex_base_with_info_2<tIndexType, K>;
using Cb_2 = CGAL::Triangulation_face_base_with_info_2<tIndexType, K>;
using Tds_2 = CGAL::Triangulation_data_structure_2<Vb_2, Cb_2>;
using DT_2 = CGAL::Delaunay_triangulation_2<K, Tds_2>;
using Point_2 = DT_2::Point;
using PointWithInfo_2 = std::pair<Point_2, tIndexType>;
using Points_2 = std::vector<PointWithInfo_2>;
using SearchTraits_2 = CGAL::Spatial_sort_traits_adapter_2<K, CGAL::First_of_pair_property_map<PointWithInfo_2>>;
//        using Fh_2 = DT_2::Face_handle;
//        using Vh_2 = DT_2::Vertex_handle;
//        using Circ_2 = CGAL::Circle_2<K>;
//        using Box_2 = CGAL::Bbox_2;
//        using Vec_2 = CGAL::Vector_2<K>;

//3D definitions
using Vb_3 = CGAL::Triangulation_vertex_base_with_info_3<tIndexType, K>;
using Cb_3 = CGAL::Triangulation_cell_base_with_info_3<tIndexType, K>;
using Tds_3 = CGAL::Triangulation_data_structure_3<Vb_3, Cb_3>;
using DT_3 = CGAL::Delaunay_triangulation_3<K, Tds_3>;
using Point_3 = DT_3::Point;
using PointWithInfo_3 = std::pair<Point_3, tIndexType>;
using Points_3 = std::vector<PointWithInfo_3>;
using SearchTraits_3 = CGAL::Spatial_sort_traits_adapter_3<K, CGAL::First_of_pair_property_map<PointWithInfo_3>>;
//        using Fh_3 = DT_3::Cell_handle;
//        using Vh_3 = DT_3::Vertex_handle;
//        using Circ_3 = CGAL::Sphere_3<K>;
//        using Box_3 = CGAL::Bbox_3;
//        using Vec_3 = CGAL::Vector_3<K>;

template<tDimType D, typename Precision, class PointArray>
std::ostream &operator<<(std::ostream &os, const Point<D, Precision, PointArray> &p) {
    os << "[" << p[0];

    for (tDimType d = 1; d < D; ++d) {
        os << ", " << p[d];

    }
    os << "]";

    return os;
}

template<tDimType D, class SimplexArray>
std::ostream &operator<<(std::ostream &os, const Simplex<D, SimplexArray> &s) {
    os << "[" << s.vertex(0);

    for (tDimType d = 1; d < D + 1; ++d) {
        os << ", " << s.vertex(d);

    }
    os << "]";
    os << " (" << s.neighbor(0);

    for (tDimType d = 1; d < D + 1; ++d) {
        os << ", " << s.neighbor(d);

    }
    os << ")";

    return os;
}


template<tDimType D, typename Precision>
struct Generator;

template<typename Precision>
struct Generator<2, Precision> {

    Points_2 generate(const tIndexType &n) {

        std::mt19937 generator(SEED);
        std::uniform_real_distribution<Precision> distribution;
        auto rand = std::bind(distribution, generator);

        Points_2 points;
        points.reserve(n);

        for (tIndexType i = 0; i < n; ++i) {
            Point_2 p_cgal(rand(), rand());
            points.push_back(std::make_pair(p_cgal, i));
        }

        SearchTraits_2 traits;
        CGAL::spatial_sort(points.begin(), points.end(), traits);

        return points;
    }

    template<class PointArray>
    void convert(PointArray &pa, const Points_2 &points) {
        pa.ensure(points.size());
        for (tIndexType i = 0; i < points.size(); ++i) {
            auto p = pa.get(i);

            p[0] = points[i].first.x();
            p[1] = points[i].first.y();
        }
    }
};

template<typename Precision>
struct Generator<3, Precision> {

    Points_3 generate(const tIndexType &n) {

        std::mt19937 generator(SEED);
        std::uniform_real_distribution<Precision> distribution;
        auto rand = std::bind(distribution, generator);

        Points_3 points;
        points.reserve(n);

        for (tIndexType i = 0; i < n; ++i) {
            Point_3 p_cgal(rand(), rand(), rand());
            points.push_back(std::make_pair(p_cgal, i));
        }

        SearchTraits_3 traits;
        CGAL::spatial_sort(points.begin(), points.end(), traits);

        return points;
    }

    template<class PointArray>
    void convert(PointArray &pa, const Points_3 &points) {
        pa.ensure(points.size());
        for (tIndexType i = 0; i < points.size(); ++i) {
            auto p = pa.get(i);

            p[0] = points[i].first.x();
            p[1] = points[i].first.y();
            p[2] = points[i].first.z();
        }
    }
};

template<tDimType D>
struct Triangulator;

template<>
struct Triangulator<2> {

    auto cgal(const Points_2 &P) {

        DT_2 T(P.begin(), P.end());

        return T;
    }

    template<class SimplexArray, class DT>
    auto convert(const DT &T) {

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

    auto cgal(const Points_3 &P) {

        DT_3 T(P.begin(), P.end());

        return T;
    }

    template<class SimplexArray, class DT>
    auto convert(const DT &T) {

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

        if constexpr (SimplexArray::isVectorized) {
            for (tIndexType i = 0;
                 i + Vc::Vector<tIndexType>::size() - 1 < simplices.size(); i += Vc::Vector<tIndexType>::size()) {

                for (tDimType d = 0; d < D + 1; ++d) {

                    auto neighbors = simplices.neighbors.vec(i, d);
                    auto mask = neighbors == INF;
                    neighbors(mask) = 0;

                    for (tDimType d2 = 0; d2 < D + 1; ++d2) {

                        auto p = simplices.vertices.vec(neighbors, d2);
                        auto det = insphere_fast<D, Precision>(i, p, simplices, points);

                        if (Vc::any_of((det < 0) & !Vc::simd_cast<Vc::Mask<Precision>>(mask))) {
                            valid = false;
                        }
                    }
                }
            }
        } else {
            for (tIndexType i = 0; i < simplices.size(); ++i) {
                auto s = simplices.get(i);

                for (tDimType n = 0; n < D + 1; ++n) {
                    if (s.neighbor(n) != INF) {
                        auto sn = simplices.get(n);

                        for (tDimType j = 0; j < D + 1; ++j) {
                            Precision det = insphere_fast<D, Precision>(i, sn.vertex(j), simplices, points);

                            if (det < 0) {
                                valid = false;
                            }
                        }
                    }
                }
            }
        }

        return valid;
    }
};

template<class SimplexArray, class PointArray>
void timeFunction(SimplexArray &simplices, const PointArray &points) {
    Checker<SimplexArray::D, typename SimplexArray::Precision> checker;

    auto t1 = std::chrono::high_resolution_clock::now();
    simplices.precompute(points);
    auto t2 = std::chrono::high_resolution_clock::now();
    auto tPrep = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

    auto t3 = std::chrono::high_resolution_clock::now();
    bool valid = checker.check(simplices, points);
    auto t4 = std::chrono::high_resolution_clock::now();
    auto tCheck = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3).count();

    std::cout << "Layout: "
              << SimplexArray::template MemoryLayout<typename SimplexArray::Precision, SimplexArray::D>::name()
              << " valid: " << valid
              << " Precomp: " << (SimplexArray::hasSubdets ? std::to_string(tPrep) : "no")
              << " Check: " << tCheck
              << " Total: " << ((SimplexArray::hasSubdets ? tPrep : 0)) + tCheck
              << std::endl;

}

template<class SimplexArray1, class PointArray1, class SimplexArray2, class PointArray2>
void verify(const SimplexArray1 &simplices1, const PointArray1 &points1,
            const SimplexArray2 &simplices2, const PointArray2 &points2) {

    std::cout << "Checking\n"
              << "\tLayout1: "
              << SimplexArray1::template MemoryLayout<typename SimplexArray1::Precision, SimplexArray1::D>::name()
              << "\tLayout2: "
              << SimplexArray2::template MemoryLayout<typename SimplexArray2::Precision, SimplexArray2::D>::name()
              << std::endl;

    // check number of points
    if (points1.size() != points2.size()) {
        std::cout << "Number points: " << points1.size() << " vs " << points2.size() << std::endl;
    }

    // check points
    for (tIndexType i = 0; i < points1.size(); ++i) {
        auto p1 = points1.get(i);
        auto p2 = points2.get(i);

        if (!(p1 == p2)) {
            std::cout << "Points " << i << " differ: (" << p1 << ") vs (" << p2 << ")" << std::endl;
        }
    }

    // check number of simplices
    if (simplices1.size() != simplices2.size()) {
        std::cout << "Number simplices: " << simplices1.size() << " vs " << simplices2.size() << std::endl;
    }

    // check simplices
    for (tIndexType i = 0; i < simplices1.size(); ++i) {
        auto s1 = simplices1.get(i);
        auto s2 = simplices2.get(i);

        if (!(s1 == s2)) {
            std::cout << "Simplices " << i << " differ: (" << s1 << ") vs (" << s2 << ")" << std::endl;
        }

        if constexpr (SimplexArray2::isVectorized) {

            // vector of indices access
            for (Vc::Vector<tIndexType> i = Vc::Vector<tIndexType>::IndexesFromZero();
                 i.max() < simplices2.size(); i += Vc::Vector<tIndexType>(Vc::Vector<tIndexType>::size())) {

                Vc::Vector<tIndexType> v[SimplexArray2::D + 1], n[SimplexArray2::D + 1];

                for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                    v[d] = simplices2.vertices.vec(i, d);
                    n[d] = simplices2.neighbors.vec(i, d);
                }

                for (std::size_t j = 0; j < Vc::Vector<tIndexType>::size(); ++j) {
                    auto s1 = simplices1.get(i[j]);
                    bool valid = true;

                    for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                        if ((s1.vertex(d) != v[d][j]) || (s1.neighbor(d) != n[d][j])) {
                            valid = false;
                        }
                    }

                    if (!valid) {
                        std::cout << "Simplices " << j << ": " << i[j]
                                  << " differ: (" << s1 << ") vs (" << "[" << v[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            std::cout << ", " << v[d][j];

                        }
                        std::cout << "]";
                        std::cout << " (" << n[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            std::cout << ", " << n[d][j];

                        }
                        std::cout << ")" << ")" << std::endl;
                    }
                }
            }

            // contiguous load access
            for (tIndexType i = 0;
                 i + Vc::Vector<tIndexType>::size() - 1 < simplices2.size(); i += Vc::Vector<tIndexType>::size()) {

                Vc::Vector<tIndexType> v[SimplexArray2::D + 1], n[SimplexArray2::D + 1];

                for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                    v[d] = simplices2.vertices.vec(i, d);
                    n[d] = simplices2.neighbors.vec(i, d);
                }

                for (std::size_t j = 0; j < Vc::Vector<tIndexType>::size(); ++j) {
                    auto s1 = simplices1.get(i + j);
                    bool valid = true;

                    for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                        if ((s1.vertex(d) != v[d][j]) || (s1.neighbor(d) != n[d][j])) {
                            valid = false;
                        }
                    }

                    if (!valid) {
                        std::cout << "Simplices " << j << ": " << i + j
                                  << " differ: (" << s1 << ") vs (" << "[" << v[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            std::cout << ", " << v[d][j];

                        }
                        std::cout << "]";
                        std::cout << " (" << n[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            std::cout << ", " << n[d][j];

                        }
                        std::cout << ")" << ")" << std::endl;
                    }
                }
            }
        }
    }

}

#define D 3
#define Precision float

#ifdef NDEBUG
#define N 1e6
#else
#define N 1e2
#endif

int main() {

#ifdef HAS_VTUNE
    __itt_pause();
#endif

#ifdef HAS_ADVISOR
    ANNOTATE_DISABLE_COLLECTION_PUSH;
#endif

    Generator<D, Precision> generator;
    auto cgal_points = generator.generate(N);

    PointArray<Traits<D, Precision, MemoryLayoutAoA, NoPrecomputation>> points_aoa;
    generator.convert(points_aoa, cgal_points);

    PointArray<Traits<D, Precision, MemoryLayoutPA, NoPrecomputation>> points_pa;
    generator.convert(points_pa, cgal_points);

#ifdef HAS_Vc
    PointArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, NoPrecomputation>> points_vaoa;
    generator.convert(points_vaoa, cgal_points);

    PointArray<Traits<D, Precision, MemoryLayoutVectorizedPA, NoPrecomputation>> points_vpa;
    generator.convert(points_vpa, cgal_points);
#endif

    Triangulator<D> triangulator;

    auto cgal_DT = triangulator.cgal(cgal_points);

    auto simplices_aoa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutAoA, NoPrecomputation>>>(
            cgal_DT);
    auto simplices_pa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutPA, NoPrecomputation>>>(
            cgal_DT);

    auto simplices_aoa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutAoA, PrecomputeSubDets>>>(
            cgal_DT);
    auto simplices_pa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutPA, PrecomputeSubDets>>>(
            cgal_DT);

#ifdef HAS_Vc
    auto simplices_vaoa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, NoPrecomputation>>>(
            cgal_DT);
    auto simplices_vpa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedPA, NoPrecomputation>>>(
            cgal_DT);
    auto simplices_vgpa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, NoPrecomputation>>>(
            cgal_DT);

    auto simplices_vaoa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, PrecomputeSubDets>>>(
            cgal_DT);
    auto simplices_vpa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedPA, PrecomputeSubDets>>>(
            cgal_DT);
    auto simplices_vgpa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, PrecomputeSubDets>>>(
            cgal_DT);
#endif

    verify(simplices_aoa_np, points_aoa, simplices_aoa_np, points_aoa);
    verify(simplices_aoa_np, points_aoa, simplices_pa_np, points_pa);

    verify(simplices_aoa_np, points_aoa, simplices_aoa_wp, points_aoa);
    verify(simplices_aoa_np, points_aoa, simplices_pa_wp, points_pa);

#ifdef HAS_Vc
    verify(simplices_aoa_np, points_aoa, simplices_vaoa_np, points_vaoa);
    verify(simplices_aoa_np, points_aoa, simplices_vpa_np, points_vpa);
    verify(simplices_aoa_np, points_aoa, simplices_vgpa_np, points_vpa);

    verify(simplices_aoa_np, points_aoa, simplices_vaoa_wp, points_vaoa);
    verify(simplices_aoa_np, points_aoa, simplices_vpa_wp, points_vpa);
    verify(simplices_aoa_np, points_aoa, simplices_vgpa_wp, points_vpa);
#endif

#ifdef HAS_VTUNE
    __itt_resume();
#endif

#ifdef HAS_ADVISOR
    ANNOTATE_DISABLE_COLLECTION_POP;
#endif

//    timeFunction(simplices_aoa_np, points_aoa);
//    timeFunction(simplices_pa_np, points_pa);
//
//    timeFunction(simplices_aoa_wp, points_aoa);
//    timeFunction(simplices_pa_wp, points_pa);

#ifdef HAS_Vc
//    timeFunction(simplices_vaoa_np, points_vaoa);
//    timeFunction(simplices_vpa_np, points_vpa);
//    timeFunction(simplices_vgpa_np, points_vpa);
//
//    timeFunction(simplices_vaoa_wp, points_vaoa);
//    timeFunction(simplices_vpa_wp, points_vpa);
//    timeFunction(simplices_vgpa_wp, points_vpa);
#endif

    return 0;

}
