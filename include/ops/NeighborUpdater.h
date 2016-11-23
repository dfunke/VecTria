#pragma once
/*
 * Update neighbors for given list of simplices + their neighboring vertices
 */

// common includes
#include <assert.h>
#include "common_types.h"

// TBBNeighborUpdater includes
#include <tbb/parallel_do.h>
#include "tbb_types.h"

struct TBBNeighborUpdater {

    template<typename DT>
    void operator()(DT &dt,
                    const tIdSet &simplices,
                    const tFacetteLookUp &facetteLookUp) {

        const auto D = DT::tSimplex::tPoint::D;

        // setup thread local storage
        tbbETS<typename DT::tConstHandle> tsDTHandle(std::ref(dt));

        // only one mark needed
        auto queuedMark = ++dt.mark;
        auto doneMark = ++dt.mark;
        assert(queuedMark < doneMark);

        tbb::parallel_do(simplices.range(), [&](const tIdType &id,
                                                tbb::parallel_do_feeder <tIdType> &feeder) {

            auto &dtHandle = tsDTHandle.local();
            auto &simplex = dtHandle[id];

            if (simplex.mark == doneMark) {
                return; //already checked
            }
            simplex.mark = doneMark;

            for (tDimType d = 0; d < D + 1; ++d) {

                /* We need to update the neighbor if
                 *  a) it is finite, but doesn't exist in the triangulation anymore or is not a neighbor anymore
                 *  b) it is infinite
                 */

                if (!typename DT::tSimplex::isFinite(simplex.neighbors[d]) // is infinite
                    ||
                    // not in triangulation anymore or not a neighbor anymore
                    (!dtHandle.contains(simplex.neighbors[d]) ||
                     !simplex.isNeighbor(dtHandle[simplex.neighbors[d]]))
                        ) {

                    // update needed
                    simplex.neighbors[d] = dSimplex<D, Precision>::cINF; // reset to infinite

                    auto facetteHash = simplex.faceFingerprint(d);
                    auto range = facetteLookUp.get(facetteHash);

                    for (auto it = range.first; it != range.second; ++it) {
                        auto cand = (*it).second;
                        if (cand != simplex.id //its not me
                            && typename DT::tSimplex::isFinite(cand) // its finite
                            && dtHandle.contains(cand) // its in the triangulation
                            && simplex.isNeighbor(dtHandle[cand]) // it is a neighbor
                            &&
                            !dtHandle[cand].contains(simplex.vertices[d]) // its neighboring at the correct edge
                                ) {

                            simplex.neighbors[d] = cand;
                            if (dtHandle[cand].mark < queuedMark) { //TODO possible race condition
                                dtHandle[cand].mark = queuedMark;
                                feeder.add(cand);
                            }
                        }
                    }
                }
            }
        });
    }
};
