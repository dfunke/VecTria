#pragma once
/*
 * Given a set of partial triangulations and an edge triangulation,
 * stitch the three together to yield one triangulation
 */

// common includes
#include <assert.h>
#include "common_types.h"
#include "ops/FacetteLookUpBuilder.h"

// TBBTriangulationStitcher includes
#include <tbb/parallel_for.h>
#include "tbb_types.h"

template<typename FacetteLookUpBuilder, typename NeighborhoodUpdater, typename Sorter>
struct TBBTriangulationStitcher {

    template<typename PartialDTs, typename DT, typename Partitioning>
    DT operator()(PartialDTs &&partialDTs,
                  DT &&edgeDT,
                  const tIdSet &edgeSimplices,
                  const Partitioning &partitioning) {

        const auto D = DT::tSimplex::tPoint::D;

        // combine the partial triangulations into one and delete all edge simplices
        auto mergedDT = std::move(partialDTs[0]);
        for (uint i = 1; i < partialDTs.size(); ++i) {
            mergedDT.merge(std::move(partialDTs[i]));
        }
        mergedDT.filter(edgeSimplices);

        // setup thread local storage
        tbbETS<typename DT::tHandle> tsMergedDTHandle(std::ref(mergedDT));
        tbbETS<typename DT::tConstHandle> tsEdgeDTHandle(std::ref(edgeDT));

        //build the facette lookup table for simplices neighboring edge simplices
        auto facetteLookUp = facetteLookUpBuilder(mergedDT, edgeSimplices);

        auto cmpSortFingerprints =
                [&tsMergedDTHandle](const tIdType &a, const tIdType &b) {
                    return tsMergedDTHandle.local()[a].fingerprint() < tsMergedDTHandle.local()[b].fingerprint();
                };

        tIdVector deletedSimplices(edgeSimplices.begin(), edgeSimplices.end());

        sorter(deletedSimplices.begin(), deletedSimplices.end(), cmpSortFingerprints);

        // setup thread local storage for facette lookup table
        tbbETS<tConcurrentGrowingFacetteLookUpHandle> tsFacetteLookUpHandle(std::ref(facetteLookUp));

        tConcurrentIdVector cConvexHull;
        tConcurrentIdVector insertedSimplices;

        auto cmpSearchFingerprints =
                [&tsMergedDTHandle](const tIdType &a, const tHashType &hash) {
                    return tsMergedDTHandle.local()[a].fingerprint() < hash;
                };

        tbb::parallel_for(edgeDT.range(), [&](auto &r) {

            auto &mergedDTHandle = tsMergedDTHandle.local();
            auto &facetteLookUpHandle = tsFacetteLookUpHandle.local();

            // we need an explicit iterator loop here for the it < r.end() comparision
            // range-based for loop uses it != r.end() which doesn't work
            for (auto it = r.begin(); it < r.end(); ++it) {
                auto &edgeSimplex = *it;

                // check whether edgeSimplex is completely contained in one partition
                uint p0 = partitioning.partition(edgeSimplex.vertices[0]);
                bool insert = !partitioning[p0].contains(edgeSimplex);

                if (!insert) {
                    // the simplex is completely in one partition -> it must have been found
                    // before
                    auto lb =
                            std::lower_bound(deletedSimplices.begin(), deletedSimplices.end(),
                                             edgeSimplex.fingerprint(), cmpSearchFingerprints);
                    assert(lb == deletedSimplices.end() || *lb != 0);

                    for (; lb != deletedSimplices.end() &&
                           mergedDTHandle[*lb].fingerprint() == edgeSimplex.fingerprint(); ++lb) {
                        if (edgeSimplex.equalVertices(mergedDTHandle[*lb])) {
                            insert = true;
                        }
                    }
                }

                if (insert) {
                    insertedSimplices.push_back(edgeSimplex.id);

                    if (!edgeSimplex.isFinite())
                        cConvexHull.push_back(edgeSimplex.id);

                    for (uint d = 0; d < D + 1; ++d) {
                        facetteLookUpHandle.insert((edgeSimplex.faceFingerprint(d)), edgeSimplex.id);
                    }
                } else {
                    edgeSimplex.id = typename DT::tSimplex::cINF;
                }
            }
        });

        mergedDT.merge(std::move(edgeDT), false);
        mergedDT.convexHull.reserve(mergedDT.convexHull.size() + cConvexHull.size());
        std::move(cConvexHull.begin(), cConvexHull.end(), std::back_inserter(mergedDT.convexHull));

        neighborhoodUpdater(mergedDT, insertedSimplices, facetteLookUp);

        return mergedDT;
    }

private:
    FacetteLookUpBuilder<FacetteLookUpPolicy::pNoIncludeSelf, FacetteLookUpPolicy::pIncludeNeighbor> facetteLookUpBuilder;
    Sorter sorter;
    NeighborhoodUpdater neighborhoodUpdater;
};
