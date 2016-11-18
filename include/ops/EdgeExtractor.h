#pragma once
/*
 * Starting from the convex hull of the DT,
 * test a simplex whether its circumsphere is within the specified polytope of the other partitions.
 * If so, enqueue all its neighbors for processing.
 */

// common includes
#include <assert.h>
#include "common_types.h"

// TBBEdgeExtractor includes
#include <tbb/parallel_do.h>
#include "tbb_types.h"

struct TBBEdgeExtractor {

    template<typename DT, typename POINTS, typename POLYTOPE>
    void getEdge(const DT &dt,
                 const POINTS &points,
                 const POLYTOPE &polytope,
                 tConcurrentGrowingIdSet &edgePoints,
                 tConcurrentGrowingIdSet &edgeSimplices) {

        // add infinite points to edge
        auto edgePointsHandle = edgePoints.handle();
        for (auto k = typename DT::tSimplex::tPoint::cINF; k != 0; ++k)
            edgePointsHandle.insert(k);

        // setup thread local storage
        tbbETS <tConcurrentGrowingIdSetHandle> tsEdgePointsHandle(std::ref(edgePoints));
        tbbETS <tConcurrentGrowingIdSetHandle> tsEdgeSimplicesHandle(std::ref(edgeSimplices));
        tbbETS<typename DT::tConstHandle> tsDTHandle(std::ref(dt));

        auto queuedMark = ++simplices.mark;
        auto doneMark = ++simplices.mark;
        assert(queuedMark < doneMark);

        tbb::parallel_do(dt.convexHull, [&](const tIdType id,
                                            tbb::parallel_do_feeder <tIdType> &feeder) {

            if (!typename DT::tSimplex::isFinite(id))
                return;

            // acquire thread local handles
            auto &edgePointsHandle = tsEdgePointsHandle.local();
            auto &edgeSimplicesHandle = tsEdgeSimplicesHandle.local();
            auto &dtHandle = tsDTHandle.local();

            const auto &simplex = dtHandle[id];

            if (simplex.mark == doneMark) {
                return; //already checked
            }
            simplex.mark = doneMark;

            const auto cs = simplex.circumsphere(points);
            if (polytope.intersects(cs)) {

                edgeSimplicesHandle.insert(simplex.id);

                for (const auto &v : simplex.vertices) {
                    edgePointsHandle.insert(v);
                }

                for (const auto &n : simplex.neighbors) {
                    if (typename DT::tSimplex::isFinite(n) && dtHandle[n].mark < queuedMark) {
                        // n was not yet inspected
                        dtHandle[n].mark = queuedMark;
                        feeder.add(n);
                    }
                }
            }
        });
    }
};
