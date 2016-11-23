#pragma once
/*
 * Build lookup table of facette hashs for given simplices + their neighboring simplices
 */

// common includes
#include <assert.h>
#include "common_types.h"

// TBBFacetteLookUpBuilder includes
#include <tbb/parallel_for.h>
#include "tbb_types.h"

template<bool INCLUDE_SELF = false, // include the simplices of idSet to lookup table
        bool INCLUDE_NEIGHBOR = true> // include the neighbors of simplices in idSet to lookup table
struct TBBFacetteLookUpBuilder {

    template<typename DT>
    tConcurrentGrowingFacetteLookUp operator()(const DT &dt,
                                               const tIdSet &simplices) {

        const auto D = DT::tSimplex::tPoint::D;

        // initialize lookup table
        tConcurrentGrowingFacetteLookUp facetteLookUp((INCLUDE_SELF + INCLUDE_NEIGHBOR) * (D + 1) * simplices.size());

        // setup thread local storage
        tbbETS <tConcurrentGrowingFacetteLookUpHandle> tsFacetteLookUpHandle(std::ref(facetteLookUp));
        tbbETS<typename DT::tConstHandle> tsDTHandle(std::ref(dt));

        // only one mark needed
        auto doneMark = ++dt.mark;

        tbb::parallel_for(simplices.range(), [&](const auto &range) {

            auto &facetteLookUpHandle = tsFacetteLookUpHandle.local();
            auto &dtHandle = tsDTHandle.local();

            for (const auto &simplexID : range) {
                const auto &simplex = dtHandle[simplexID];

                if (simplex.mark == doneMark) {
                    continue; //already checked
                }
                simplex.mark = doneMark;

                if (INCLUDE_SELF) {
                    for (tDimType i = 0; i < D + 1; ++i) {
                        auto facetteHash = simplex.faceFingerprint(i);
                        facetteLookUpHandle.insert(facetteHash, simplexID);
                    }
                }

                if (INCLUDE_NEIGHBOR) {
                    for (const auto &n : simplex.neighbors) {
                        if (typename DT::tSimplex::isFinite(n) && !simplices.count(n)) {
                            // we have an "inward" neighbor, add it to the lookup table
                            if (dtHandle[n].mark < doneMark) { //TODO possible race condition
                                dtHandle[n].mark = doneMark;

                                for (tDimType i = 0; i < D + 1; ++i) {
                                    auto facetteHash = dtHandle[n].faceFingerprint(i);
                                    facetteLookUpHandle.insert(facetteHash, n);
                                }
                            }
                        }
                    }
                }
            }
        });

        return facetteLookUp;
    }
};
