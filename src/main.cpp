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

#include "Stats.h"

#include "Generator.h"
#include "Triangulator.h"
#include "Checker.h"

#include "FileReader.h"
#include "TextTable.h"

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
        STAT_CLEAR;
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
    output.add(std::to_string(STAT_GET(Insphere)));
    output.add(std::to_string(STAT_GET(StaticFilterFail)));
    output.add(std::to_string(STAT_GET(PermanentFill)));
    output.add(std::to_string(STAT_GET(PermanentFilterFail)));
    output.add(std::to_string(STAT_GET(AdaptiveFilterFail)));
    output.endOfRow();

}

template<class PointArray>
void timeDeltaCalc(TextTable &output) {

    FileReader<PointArray::D, typename PointArray::Precision> freader;
    PointArray p1, p2;

    // const std::string DIR = "/home/funke/devel/KineticDelaunay/data/crack/";
    const std::string DIR = "/Users/funke/Development/ParDeTria/data/crack/";

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

template<class PointArray, class SimplexArray>
void timeKineticValid(TextTable &output) {

    output.add("Layout");
    output.add("AvgDelta");
    output.add("#points");
    output.add("#simplices");
    output.add("% invalid");
    output.add("Time Check");
    output.endOfRow();

    FileReader<PointArray::D, typename PointArray::Precision> freader;
    PointArray pa, pb;
    Points_3 p_cgal;

    Triangulator<PointArray::D> triangulator;
    Generator<PointArray::D, typename PointArray::Precision> conv;
    Checker<SimplexArray::D, typename SimplexArray::Precision> checker;

    // const std::string DIR = "/home/funke/devel/KineticDelaunay/data/crack/";
    const std::string DIR = "/Users/funke/Development/ParDeTria/data/crack/";

    std::filesystem::path p(DIR);
    std::filesystem::directory_iterator start(p);
    std::filesystem::directory_iterator end;

    std::vector<std::filesystem::path> files;
    std::copy(start, end, std::back_inserter(files));
    std::sort(files.begin(), files.end());

    freader.read(files[0], pa);
    conv.convert(p_cgal, pa);

    auto dt_cgal = triangulator.cgal(p_cgal);
    auto dt = triangulator.template convert<SimplexArray>(dt_cgal);

//    for (uint i = 1; i < files.size(); ++i) {
    for (uint i = 1; i < 6; ++i) {
        std::cout << "Processing file " << files[i] << std::endl;
        freader.read(files[i], pb);

        auto delta = pa.getDelta(pb);

        double avg = 0;
        for (tIndexType i = 0; i < delta.size(); ++i) {
            for (tDimType d = 0; d < PointArray::D; ++d) {
                avg += delta(i, d);
            }
        }
        avg /= (PointArray::D * delta.size());

        auto t1 = std::chrono::high_resolution_clock::now();
        auto invalid = checker.check(dt, pb);
        auto t2 = std::chrono::high_resolution_clock::now();

        conv.convert(p_cgal, pb);
        pa = pb;
        dt_cgal = triangulator.cgal(p_cgal);
        dt = triangulator.template convert<SimplexArray>(dt_cgal);

        output.add(PointArray::template MemoryLayout<typename PointArray::Precision, PointArray::D>::name());
        output.add(std::to_string(avg));
        output.add(std::to_string(pa.size()));
        output.add(std::to_string(dt.size()));
        output.add(std::to_string(invalid / static_cast<double>(dt.size())));
        output.add(std::to_string(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()));
        output.endOfRow();
    }

}

template<class PointArray, class SimplexArray, class Movement, typename... Args>
void timeSyntheticKineticValid(TextTable &output,
                               const tIndexType N, const tIndexType T,
                               Args... args) {
 
    output.add("Mover");
    output.add("AvgDelta");
    output.add("#points");
    output.add("#simplices");
    output.add("% invalid");
    output.add("Time Check");
    output.endOfRow();

    PointArray pa, pb;
    Points_3 p_cgal;

    Triangulator<PointArray::D> triangulator;
    Generator<PointArray::D, typename PointArray::Precision> generator;
    Checker<SimplexArray::D, typename SimplexArray::Precision> checker;
    Movement mover;

    p_cgal = generator.generate(N);
    generator.convert(pa, p_cgal);
    pb = pa;
    Predicates<typename PointArray::Precision>::set_static_limits(1, 1, 1);

    auto dt_cgal = triangulator.cgal(p_cgal);
    auto dt = triangulator.template convert<SimplexArray>(dt_cgal);
    
    auto invalid = checker.check(dt, pa);
    assert(!invalid);

    for (uint i = 1; i < T; ++i) {
        mover.move(pb, args...);

        auto delta = pa.getDelta(pb);

        double avg = 0;
        for (tIndexType i = 0; i < delta.size(); ++i) {
            for (tDimType d = 0; d < PointArray::D; ++d) {
                avg += delta(i, d);
            }
        }
        avg /= (PointArray::D * delta.size());

        auto t1 = std::chrono::high_resolution_clock::now();
        auto invalid = checker.check(dt, pb);
        auto t2 = std::chrono::high_resolution_clock::now();

        generator.convert(p_cgal, pb);
        pa = pb;
        dt_cgal = triangulator.cgal(p_cgal);
        dt = triangulator.template convert<SimplexArray>(dt_cgal);

        output.add(Movement::name());
        output.add(std::to_string(avg));
        output.add(std::to_string(pa.size()));
        output.add(std::to_string(dt.size()));
        output.add(std::to_string(invalid / static_cast<double>(dt.size())));
        output.add(std::to_string(std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count()));
        output.endOfRow();
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

    using tTraits = Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, NoOppVertex>;
    
    Precision stepSize = 0.01;
    tVector<D, Precision> step = {{stepSize, stepSize, stepSize}};
    
    TextTable output;
    
    timeSyntheticKineticValid<PointArray<tTraits>, SimplexArray<tTraits>, UniformMover>(output, N, 10, step);
    
    output.endOfRow();
    
    timeSyntheticKineticValid<PointArray<tTraits>, SimplexArray<tTraits>, UniformStochasticMover>(output, N, 10, step);
    
    std::vector<Box<tVector<D, Precision>>> zones;
    std::vector<tVector<D, Precision>> steps;
    
    zones.push_back(Box<tVector<D, Precision>>({{0.0, 0.0, 0.0}}, {{0.5, 0.5, 1.0}}));
    steps.push_back({{stepSize, -stepSize, 0.0}});
    
    zones.push_back(Box<tVector<D, Precision>>({{0.5, 0.0, 0.0}}, {{1.0, 0.5, 1.0}}));
    steps.push_back({{stepSize, stepSize, 0.0}});
    
    zones.push_back(Box<tVector<D, Precision>>({{0.5, 0.5, 0.0}}, {{1.0, 1.0, 1.0}}));
    steps.push_back({{-stepSize, stepSize, 0.0}});
    
    zones.push_back(Box<tVector<D, Precision>>({{0.0, 0.5, 0.0}}, {{0.5, 1.0, 1.0}}));
    steps.push_back({{-stepSize, -stepSize, 0.0}});
    
    output.endOfRow();
    
    timeSyntheticKineticValid<PointArray<tTraits>, SimplexArray<tTraits>, ZonedMover>(output, N, 10, zones, steps);
    
    output.endOfRow();
    
    //timeSyntheticKineticValid<PointArray<tTraits>, SimplexArray<tTraits>, ZonedStochasticMover>(output, N, 10, zones, steps);
    
    std::cout << output;

    
    return 0;

}
