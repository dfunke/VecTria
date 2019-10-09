#pragma once

#include "Datastructures.h"
#include "GeometryStructures.h"
#include "CGAL.h"

template<tDimType D, typename Precision>
struct Checker;

template<typename Precision>
struct Checker<3, Precision> {

    static constexpr tDimType D = 3;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// vectorized

    template<class SimplexArray, class PointArray, typename Ret=bool>
    auto check(const SimplexArray &simplices,
               const PointArray &points) -> std::enable_if_t<SimplexArray::isVectorized, Ret> {

        bool valid = true;


        for (tIndexType i = 0;
             i + Vc::Vector<tIndexType>::size() - 1 < simplices.size(); i += Vc::Vector<tIndexType>::size()) {

            for (tDimType d = 0; d < D + 1; ++d) {

                Vc::Vector<tIndexType> oppVertex(INF);

                if constexpr (SimplexArray::hasOppVertex) {
                    oppVertex = simplices.opp_vertex.vec(i, d);
                } else {

                    auto s = Vc::Vector<tIndexType>(i) + Vc::Vector<tIndexType>::IndexesFromZero();

                    auto neighbors = simplices.neighbors.vec(i, d);
                    auto maskInfNeighbor = neighbors == INF;
                    neighbors(maskInfNeighbor) = 0;

                    for (tDimType d1 = 0; d1 < D + 1; ++d1) {

                        auto n = simplices.neighbors.vec(neighbors, d1);
                        auto maskEqualNeighbor = n == s;

                        oppVertex(maskEqualNeighbor) = simplices.vertices.vec(neighbors, d1);
                    }
                }

                auto maskInfOppVertex = oppVertex == INF;
                oppVertex(maskInfOppVertex) = 0;

                auto det = Predicates<Precision>::insphere_fast(i, oppVertex, simplices, points);

                static_assert(Vc::Vector<tIndexType>::size() == Vc::Vector<Precision>::size());
                if (Vc::any_of((det < 0) & !Vc::simd_cast<Vc::Mask<Precision>>(maskInfOppVertex))) {
                    valid = false;
                }
            }
        }

        return valid;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// non-vectorized

    template<class SimplexArray, class PointArray, typename Ret=bool>
    auto check(const SimplexArray &simplices,
               const PointArray &points) -> std::enable_if_t<not SimplexArray::isVectorized, Ret> {

        bool valid = true;

        for (tIndexType i = 0; i < simplices.size(); ++i) {
            for (tDimType d = 0; d < D + 1; ++d) {

                tIndexType oppVertex = INF;

                if constexpr (SimplexArray::hasOppVertex) {
                    oppVertex = simplices.opp_vertex(i, d);
                } else {
                    auto s = simplices.get(i);

                    if (s.neighbor(d) != INF) {
                        auto sn = simplices.get(d);

                        for (tDimType d1 = 0; d1 < D + 1; ++d1) {
                            if (sn.neighbor(d1) == i) {
                                oppVertex = sn.vertex(d1);
                                break;
                            }
                        }
                    }
                }

                if (oppVertex != INF) {
                    auto[det, permanent] = Predicates<Precision>::insphere(i, oppVertex, simplices, points);

                    if (!Predicates<Precision>::isFinal(det, permanent)) {
                        det = Predicates<Precision>::insphere_adapt(permanent, i, oppVertex, simplices, points);
                    }

                    if (det < 0) {
                        valid = false;
                    }
                }
            }
        }

        return valid;
    }
};
