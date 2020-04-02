#pragma once

#include "GeometryStructures.h"
#include "Triangulator.h"
#include "Generator.h"
#include "Predicates.h"

template<typename PointArray, typename SimplexArray>
struct SimplexOrientationFixer {
    using Precision = typename PointArray::Precision;
    using Pred = Predicates<Precision>;

    SimplexOrientationFixer(const PointArray *points, SimplexArray *triangulation)
            : points(points),
              triangulation(triangulation) {
    }

    void operator()() {
        for (tIndexType i = 0; i < triangulation->size(); ++i) {
            fixSimplex(i);
        }
    }

    void fixSimplex(tIndexType simplexId) const {
        if (calculateSimplexOrientation(simplexId) < 0) {
            invertSimplex(simplexId);
        }
    }

    Precision calculateSimplexOrientation(tIndexType simplexId) const {
        return Pred::pred_orient(simplexId, *triangulation, *points);
    }

    void invertSimplex(tIndexType simplexId) const {
        auto simplex = triangulation->get(simplexId);
        std::swap(simplex.vertex(0), simplex.vertex(1));
        std::swap(simplex.neighbor(0), simplex.neighbor(1));
    }

private:
    const PointArray *points;
    SimplexArray *triangulation;
};