#pragma once
/*
 * Given a set of partial triangulations and an edge triangulation,
 * stitch the three together to yield one triangulation
 */

// common includes
#include <assert.h>
#include "common_types.h"
#include "TriangulationMerger.h"

// TBBDACTriangulator includes
#include <tbb/parallel_for.h>
#include "tbb_types.h"

struct DACTriangulatorConfig {

    tRecDepthType cMaxRecursionDepth;
    tIdType cBasecaseCutoff;

};

template<typename Partitioner, typename Merger, typename BasecaseTriangulator>
struct TBBDACTriangulator {

    template<typename PointRange, typename Polytope, typename Points, typename DT>
    auto operator()(const PointRange &partitionPoints,
                    const Polytope &partitionBounds,
                    const Points &points,
                    const tRecDepthType &recDepth) {


        if (recDepth >= conf.cMaxRecursionDepth || partitionPoints.size() < conf.cBasecaseCutoff) {
            // recursion base case
            return basecaseTriangulator(partitionPoints, partitionBounds, points);
        }

        const auto partitioning = partitioner(partitionPoints, partitionBounds, points);

        std::array<DT, Partitioner::cNumPartitions> partialDTs;
        tbb::parallel_for(
                tDimType(0), Partitioner::cNumPartitons,
                [&](const tDimType i) {
                    partialDTs[i] = this->operator()(partitioning[i].points, partitioning[i].bounds,
                                                     points, recDepth + 1);
                }
        );

        return merger(std::move(partialDTs), partitioning);

    }

private:
    Partitioner partitioner;
    Merger merger;
    BasecaseTriangulator basecaseTriangulator;

    DACTriangulatorConfig conf;
};
