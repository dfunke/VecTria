#pragma once

#include <gtest/gtest.h>
#include "Generator.h"
#include "Triangulator.h"
#include "Predicates.h"
#include "GeometryStructures.h"

struct TriangulatorTestContext {
    static constexpr tDimType D = 3;
    using Precision = double;
    using DefaultTraits = Traits<D, Precision, MemoryLayoutPA, NoPrecomputation, NoOppVertex>;
    using PointArrayT = PointArray<DefaultTraits>;
    using SimplexArrayT = SimplexArray<DefaultTraits>;
    using Pred = Predicates<Precision>;

    PointArrayT points;
    SimplexArrayT delaunayTriangulation;
    tIndexType simplexIndex;
    tIndexType neighborSimplexIndex;
    tIndexType farVertexPosition;

    TriangulatorTestContext() {
        Generator<3, double> pointGenerator;
        auto cgalPoints = pointGenerator.generate(100);

        Triangulator<3> triangulator;
        auto cgalTriangulation = triangulator.cgal(cgalPoints);

        pointGenerator.convert(points, cgalPoints);
        delaunayTriangulation = triangulator.convert<SimplexArrayT>(cgalTriangulation);

        simplexIndex = 25;
        tDimType neighborPosition = 0;
        assert(simplexIndex < delaunayTriangulation.size());
        neighborSimplexIndex = delaunayTriangulation.get(simplexIndex).neighbor(neighborPosition);
        assert(neighborPosition != INF);
        farVertexPosition = delaunayTriangulation.get(simplexIndex).vertex(neighborPosition);
    }
};

TEST(Triangulator, SimplexOrientation) {
    TriangulatorTestContext context;

    auto orientation = TriangulatorTestContext::Pred::pred_orient(context.simplexIndex,
                                                                  context.delaunayTriangulation,
                                                                  context.points);
    EXPECT_GT(orientation, 0);
}


TEST(Triangulator, SimplexDelaunayProperty) {
    TriangulatorTestContext context;

    auto adjacentSimplexId = context.neighborSimplexIndex;
    auto farVertexId = context.farVertexPosition;
    auto insphereResult = TriangulatorTestContext::Pred::pred_insphere(adjacentSimplexId, farVertexId,
                                                                       context.delaunayTriangulation,
                                                                       context.points);
    EXPECT_LE(insphereResult, 0);
}