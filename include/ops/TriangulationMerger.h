#pragma once
/*
 * Given a set of partial triangulations and an edge triangulation,
 * stitch the three together to yield one triangulation
 */

// common includes
#include <assert.h>
#include "common_types.h"

// TBBTriangulationMerger includes
#include <tbb/parallel_for.h>
#include "tbb_types.h"

template<typename Triangulator, typename EdgeExtractor, typename TriangulationStitcher>
struct TBBTriangulationMerger {

    template<typename PartialDTs, typename Partitioning>
    auto operator()(PartialDTs &&partialDTs,
                    const Partitioning &partitioning) {


        tConcurrentGrowingIdSet edgePointIds/*TODO (size)*/;
        tConcurrentGrowingIdSet edgeSimplexIds/*TODO (size)*/;

        tbb::parallel_for(
                typename PartialDTs::size_type(0), partitioning.size(),
                [&](const typename PartialDTs::size_type i) {
                    edgeExtractor(partialDTs[i],
                                  partitioning.bounds - partitioning[i].bounds,
                                  edgePointIds, edgeSimplexIds);
                }
        );

        auto edgeDT = triangulator(tIdSet(std::move(edgePointIds)), partitioning.bounds);

        return triangulationStitcher(std::move(partialDTs),
                                     std::move(edgeDT),
                                     tIdSet(std::move(edgeSimplexIds)),
                                     partitioning);
    }

private:
    Triangulator triangulator;
    EdgeExtractor edgeExtractor;
    TriangulationStitcher triangulationStitcher;
};
