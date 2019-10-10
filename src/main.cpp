//
// Created by funke on 11/23/16.
//

#include <functional>
#include <chrono>

#ifdef HAS_VTUNE
#include <ittnotify.h>
#endif

#ifdef HAS_ADVISOR
#include <advisor-annotate.h>
#endif

#include "Generator.h"
#include "Triangulator.h"
#include "Checker.h"


#include "FileReader.h"
#include "TextTable.h"

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

template<class SimplexArray, class PointArray>
void
timeValidityCheck(SimplexArray &simplices, const PointArray &points, const unsigned short reps, TextTable &output) {
    Checker<SimplexArray::D, typename SimplexArray::Precision> checker;

    auto t1 = std::chrono::high_resolution_clock::now();
    for (uint i = 0; i < reps; ++i) {
        simplices.precompute(points);
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    auto tPrep = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count() / reps;

    auto t3 = std::chrono::high_resolution_clock::now();
    bool valid = false;

    for (uint i = 0; i < reps; ++i) {
        valid = checker.check(simplices, points);
    }

    auto t4 = std::chrono::high_resolution_clock::now();
    auto tCheck = std::chrono::duration_cast<std::chrono::duration<double>>(t4 - t3).count() / reps;

    output.add(SimplexArray::template MemoryLayout<typename SimplexArray::Precision, SimplexArray::D>::name());
    output.add(SimplexArray::hasOppVertex ? "yes" : "no");
    output.add(std::to_string(valid));
    output.add(SimplexArray::hasSubdets ? std::to_string(tPrep) : "no");
    output.add(std::to_string(tCheck));
    output.add(std::to_string((SimplexArray::hasSubdets ? tPrep : 0) + tCheck));
    output.endOfRow();

}

template<class PointArray>
void timeDeltaCalc(TextTable &output) {

    FileReader<PointArray::D, typename PointArray::Precision> freader;
    PointArray p1, p2;

    const std::string DIR = "/home/funke/devel/KineticDelaunay/data/crack/";

    std::filesystem::path p(DIR);
    std::filesystem::directory_iterator start(p);
    std::filesystem::directory_iterator end;

    std::vector<std::filesystem::path> files;
    std::copy(start, end, std::back_inserter(files));
    std::sort(files.begin(), files.end());

    freader.read(files[0], p1);

    double time = 0;
    double avgDelta = 0;

    for (uint i = 1; i < files.size(); ++i) {
        freader.read(files[i], p2);

        auto t1 = std::chrono::high_resolution_clock::now();
        auto delta = p1.getDelta(p2);
        auto t2 = std::chrono::high_resolution_clock::now();

        time += std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

        double avg = 0;
        for (tIndexType i = 0; i < delta.size(); ++i) {
            for (tDimType d = 0; d < PointArray::D; ++d) {
                avg += delta(i, d);
            }
        }
        avgDelta += avg / (PointArray::D * delta.size());

        p1 = p2;
    }

    output.add(PointArray::template MemoryLayout<typename PointArray::Precision, PointArray::D>::name());
    output.add(std::to_string(time));
    output.add(std::to_string(avgDelta));
    output.endOfRow();

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

#define D 3
#define Precision float

#ifdef NDEBUG
#define N 1e6
#define REPS 1
#else
#define N 1e2
#define REPS 1
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
    Predicates<Precision>::set_static_limits(1,1,1);

    PointArray<Traits<D, Precision, MemoryLayoutAoA, NoPrecomputation, NoOppVertex>> points_aoa;
    generator.convert(points_aoa, cgal_points);

    PointArray<Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, NoOppVertex>> points_pa;
    generator.convert(points_pa, cgal_points);

#ifdef HAS_Vc
    PointArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, NoPrecomputation, NoOppVertex>> points_vaoa;
    generator.convert(points_vaoa, cgal_points);

    PointArray<Traits<D, Precision, MemoryLayoutVectorizedPA, NoPrecomputation, NoOppVertex>> points_vpa;
    generator.convert(points_vpa, cgal_points);
#endif

    Triangulator<D> triangulator;

    auto t1 = std::chrono::high_resolution_clock::now();
    auto cgal_DT = triangulator.cgal(cgal_points);
    auto t2 = std::chrono::high_resolution_clock::now();
    std::cout << "Triangulation time: " << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()
              << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
    auto b = cgal_DT.is_valid();
    t2 = std::chrono::high_resolution_clock::now();
    std::cout << "CGAL verification time: "
              << std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()
              << " is " << (b ? "" : "NOT ") << "valid" << std::endl;


    auto simplices_aoa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutAoA, NoPrecomputation, NoOppVertex>>>(
            cgal_DT);
    auto simplices_pa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, NoOppVertex>>>(
            cgal_DT);

    auto simplices_aoa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutAoA, PrecomputeSubDets, NoOppVertex>>>(
            cgal_DT);
    auto simplices_pa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutPA, PrecomputeSubDets, NoOppVertex>>>(
            cgal_DT);

    auto simplices_aoa_np_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutAoA, NoPrecomputation, WithOppVertex>>>(
            cgal_DT);
    auto simplices_pa_np_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, WithOppVertex>>>(
            cgal_DT);

    auto simplices_aoa_wp_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutAoA, PrecomputeSubDets, WithOppVertex>>>(
            cgal_DT);
    auto simplices_pa_wp_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutPA, PrecomputeSubDets, WithOppVertex>>>(
            cgal_DT);

#ifdef HAS_Vc
    auto simplices_vaoa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, NoPrecomputation, NoOppVertex>>>(
            cgal_DT);
    auto simplices_vpa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedPA, NoPrecomputation, NoOppVertex>>>(
            cgal_DT);
    auto simplices_vgpa_np = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, NoPrecomputation, NoOppVertex>>>(
            cgal_DT);

    auto simplices_vaoa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, PrecomputeSubDets, NoOppVertex>>>(
            cgal_DT);
    auto simplices_vpa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedPA, PrecomputeSubDets, NoOppVertex>>>(
            cgal_DT);
    auto simplices_vgpa_wp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, PrecomputeSubDets, NoOppVertex>>>(
            cgal_DT);

    auto simplices_vaoa_np_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, NoPrecomputation, WithOppVertex>>>(
            cgal_DT);
    auto simplices_vpa_np_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedPA, NoPrecomputation, WithOppVertex>>>(
            cgal_DT);
    auto simplices_vgpa_np_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, NoPrecomputation, WithOppVertex>>>(
            cgal_DT);

    auto simplices_vaoa_wp_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, PrecomputeSubDets, WithOppVertex>>>(
            cgal_DT);
    auto simplices_vpa_wp_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedPA, PrecomputeSubDets, WithOppVertex>>>(
            cgal_DT);
    auto simplices_vgpa_wp_opp = triangulator.convert<SimplexArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, PrecomputeSubDets, WithOppVertex>>>(
            cgal_DT);

#endif

#ifdef HAS_VTUNE
    __itt_resume();
#endif

#ifdef HAS_ADVISOR
    ANNOTATE_DISABLE_COLLECTION_POP;
#endif

    TextTable output;
    output.add("Layout");
    output.add("OppVertex");
    output.add("valid");
    output.add("Precomp");
    output.add("Check");
    output.add("Total");
    output.endOfRow();

    timeValidityCheck(simplices_aoa_np, points_aoa, REPS, output);
    timeValidityCheck(simplices_pa_np, points_pa, REPS, output);

    timeValidityCheck(simplices_aoa_wp, points_aoa, REPS, output);
    timeValidityCheck(simplices_pa_wp, points_pa, REPS, output);

    timeValidityCheck(simplices_aoa_np_opp, points_aoa, REPS, output);
    timeValidityCheck(simplices_pa_np_opp, points_pa, REPS, output);

    timeValidityCheck(simplices_aoa_wp_opp, points_aoa, REPS, output);
    timeValidityCheck(simplices_pa_wp_opp, points_pa, REPS, output);

#ifdef HAS_Vc
    timeValidityCheck(simplices_vaoa_np, points_vaoa, REPS, output);
    timeValidityCheck(simplices_vpa_np, points_vpa, REPS, output);
    timeValidityCheck(simplices_vgpa_np, points_vpa, REPS, output);

    timeValidityCheck(simplices_vaoa_wp, points_vaoa, REPS, output);
    timeValidityCheck(simplices_vpa_wp, points_vpa, REPS, output);
    timeValidityCheck(simplices_vgpa_wp, points_vpa, REPS, output);

    timeValidityCheck(simplices_vaoa_np_opp, points_vaoa, REPS, output);
    timeValidityCheck(simplices_vpa_np_opp, points_vpa, REPS, output);
    timeValidityCheck(simplices_vgpa_np_opp, points_vpa, REPS, output);

    timeValidityCheck(simplices_vaoa_wp_opp, points_vaoa, REPS, output);
    timeValidityCheck(simplices_vpa_wp_opp, points_vpa, REPS, output);
    timeValidityCheck(simplices_vgpa_wp_opp, points_vpa, REPS, output);
#endif

    std::cout << output;


//    output.clear();
//    output.add("Layout");
//    output.add("Time");
//    output.add("AvgDelta");
//    output.endOfRow();
//
//    timeDeltaCalc<PointArray<Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, NoOppVertex>>>(output);
//    timeDeltaCalc<PointArray<Traits<D, Precision, MemoryLayoutAoA, NoPrecomputation, NoOppVertex>>>(output);
//
//#ifdef HAS_Vc
//    timeDeltaCalc<PointArray<Traits<D, Precision, MemoryLayoutVectorizedPA, NoPrecomputation, NoOppVertex>>>(output);
//    timeDeltaCalc<PointArray<Traits<D, Precision, MemoryLayoutVectorizedGroupedPA, NoPrecomputation, NoOppVertex>>>(
//            output);
//    timeDeltaCalc<PointArray<Traits<D, Precision, MemoryLayoutVectorizedAoA, NoPrecomputation, NoOppVertex>>>(output);
//#endif
//
//    std::cout << output;

    return 0;

}
