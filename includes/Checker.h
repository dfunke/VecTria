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

    template<class SimplexArray, class PointArray, typename Ret=tIndexType>
    auto check(const SimplexArray &simplices,
               const PointArray &points) -> std::enable_if_t<SimplexArray::isVectorized, Ret> {

        tIndexType invalid = 0;

        tIndexType i = 0;
        for (;
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

                auto[det, permanent] = Predicates<Precision>::insphere(i, oppVertex, simplices, points);

                static_assert(Vc::Vector<tIndexType>::size() == Vc::Vector<Precision>::size());

                auto maskNonFinal = (permanent != 0) & !Vc::simd_cast<Vc::Mask<Precision>>(maskInfOppVertex);
                if (Vc::any_of(maskNonFinal)) {
                    for (uint v = 0; v < Vc::Mask<Precision>::size(); ++v) {
                        if (maskNonFinal[v]) {
                            det[v] = Predicates<Precision>::insphere_adapt(permanent[v], i + v, oppVertex[v], simplices,
                                                                           points);
                        }
                    }
                }

                if (Vc::any_of((det > 0) & !Vc::simd_cast<Vc::Mask<Precision>>(maskInfOppVertex))) {
                    invalid += ((det > 0) & !Vc::simd_cast<Vc::Mask<Precision>>(maskInfOppVertex)).count();
                }

            }
        }

        for (; i < simplices.size(); ++i) {
            if (!check_simplex(i, simplices, points)) {
                ++invalid;
            }
        }

        return invalid;
    }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// non-vectorized

    template<class SimplexArray, class PointArray, typename Ret=std::pair<tIndexType, tIndexType>>
    auto check(const SimplexArray &simplices,
               const PointArray &points) -> std::enable_if_t<not SimplexArray::isVectorized, Ret> {

        tIndexType invalid = 0;
        tIndexType wrongOrient = 0;

        for (tIndexType i = 0; i < simplices.size(); ++i) {

            auto[valid, orient] = check_simplex(i, simplices, points);
            invalid += !valid;
            wrongOrient += !orient;
        }

        return std::make_pair(invalid, wrongOrient);
    }

    template<class SimplexArray, class PointArray>
    inline std::pair<bool, bool>
    check_simplex(const tIndexType &i, const SimplexArray &simplices, const PointArray &points) {

        bool bOrient;

        for (tDimType d = 0; d < D + 1; ++d) {
            tIndexType oppVertex = INF;

            if constexpr (SimplexArray::hasOppVertex) {
                oppVertex = simplices.opp_vertex(i, d);
            } else {
                auto s = simplices.get(i);

                if (s.neighbor(d) != INF) {
                    auto sn = simplices.get(s.neighbor(d));

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
                auto orient = Predicates<Precision>::pred_orient(i, simplices, points);
                bOrient = (orient > 0);

                if (permanent) {
                    det = Predicates<Precision>::insphere_adapt(permanent, i, oppVertex, simplices, points);
                }

                if (orient * det < 0) {
                    return std::make_pair(false, bOrient);
                }
            }
        }

        return std::make_pair(true, bOrient);

    }

};
