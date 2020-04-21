#pragma once

#include <algorithm>
#include <gtest/gtest.h>
#include "Generator.h"
#include "Triangulator.h"
#include "Predicates.h"
#include "GeometryStructures.h"
#include "SimplexOrientationFixer.h"

struct ComprehensiveTriangulatorOrientationTestContext {
    static constexpr tDimType D = 3;
    using Precision = double;
    using DefaultTraits = Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, NoOppVertex>;
    using PointArrayT = PointArray<DefaultTraits>;
    using SimplexArrayT = SimplexArray<DefaultTraits>;
    using Pred = Predicates<Precision>;

    PointArrayT points;
    SimplexArrayT delaunayTriangulation;

    ComprehensiveTriangulatorOrientationTestContext() {
        Generator<3, double> pointGenerator;
        auto cgalPoints = pointGenerator.generate(100);

        Triangulator<3> triangulator;
        auto cgalTriangulation = triangulator.cgal(cgalPoints);

        pointGenerator.convert(points, cgalPoints);
        delaunayTriangulation = triangulator.convert<SimplexArrayT>(cgalTriangulation);
    }

    std::vector<tIndexType> findNegativelyOrientedSimplices() const {
        std::vector<tIndexType> negativelyOrientedSimplices;
        for (tIndexType id = 0; id < delaunayTriangulation.size(); ++id) {
            auto orientation = TriangulatorTestContext::Pred::pred_orient(id,
                                                                          delaunayTriangulation,
                                                                          points);
            if (orientation < 0) {
                negativelyOrientedSimplices.push_back(id);
            }
        }

        return negativelyOrientedSimplices;
    }
};

TEST(ComprehensiveTriangulatorOrientationTest, SimplexOrientations) {
    ComprehensiveTriangulatorOrientationTestContext context;

    EXPECT_EQ(context.findNegativelyOrientedSimplices(), std::vector<tIndexType>{});
}