#pragma once

#include "Datastructures.h"
#include "GeometryStructures.h"
#include "CGAL.h"

template<tDimType D, typename Precision>
struct Checker;

template<typename Precision>
struct Checker<3, Precision> {

    static constexpr tDimType D = 3;

    template<class SimplexArray, class PointArray>
    bool check(const SimplexArray &simplices, const PointArray &points) {

        bool valid = true;

        if constexpr (SimplexArray::isVectorized) {

            if constexpr (SimplexArray::hasOppVertex) {
                for (tIndexType i = 0;
                     i + Vc::Vector<tIndexType>::size() - 1 < simplices.size(); i += Vc::Vector<tIndexType>::size()) {
                    for (tDimType d = 0; d < D + 1; ++d) {

                        auto p = simplices.opp_vertex.vec(i, d);
                        auto mask = p == INF;
                        p(mask) = 0;

                        auto det = Predicates<Precision>::insphere_fast(i, p, simplices, points);

                        static_assert(Vc::Vector<tIndexType>::size() == Vc::Vector<Precision>::size());
                        if (Vc::any_of((det < 0) & !Vc::simd_cast<Vc::Mask<Precision>>(mask))) {
                            valid = false;
                        }
                    }
                }
            } else {

                for (tIndexType i = 0;
                     i + Vc::Vector<tIndexType>::size() - 1 < simplices.size(); i += Vc::Vector<tIndexType>::size()) {

                    for (tDimType d = 0; d < D + 1; ++d) {

                        auto neighbors = simplices.neighbors.vec(i, d);
                        auto mask = neighbors == INF;
                        neighbors(mask) = 0;

                        for (tDimType d2 = 0; d2 < D + 1; ++d2) {

                            auto p = simplices.vertices.vec(neighbors, d2);
                            auto det = Predicates<Precision>::insphere_fast(i, p, simplices, points);

                            if (Vc::any_of((det < 0) & !Vc::simd_cast<Vc::Mask<Precision>>(mask))) {
                                valid = false;
                            }
                        }
                    }
                }
            }

        } else {

            if constexpr (SimplexArray::hasOppVertex) {
                for (tIndexType i = 0; i < simplices.size(); ++i) {
                    for (tDimType d = 0; d < D + 1; ++d) {

                        auto p = simplices.opp_vertex(i, d);

                        if (p != INF) {
                            auto[det, permanent] = Predicates<Precision>::insphere(i, p, simplices, points);

                            if (!Predicates<Precision>::isFinal(det, permanent)) {
                                std::cout << "adapt called" << std::endl;
                                det = Predicates<Precision>::insphere_adapt(permanent, i, p, simplices, points);
                            }

                            det = Predicates<Precision>::pred_orient(i, simplices, points) *
                                  Predicates<Precision>::pred_insphere(i, p, simplices, points);

                            if (det >= 0) {
//                                std::cout << "invalid" << std::endl;
                                valid = false;
                            }
                        }
                    }
                }
            } else {
                uint inv = 0;
                uint val = 0;
                for (tIndexType i = 0; i < simplices.size(); ++i) {
                    auto s = simplices.get(i);

                    for (tDimType n = 0; n < D + 1; ++n) {
                        if (s.neighbor(n) != INF) {
                            auto sn = simplices.get(n);

                            for (tDimType j = 0; j < D + 1; ++j) {
                                auto[det, permanent] = Predicates<Precision>::insphere(i, sn.vertex(j), simplices,
                                                                                       points);

                                if (!Predicates<Precision>::isFinal(det, permanent)) {
                                    std::cout << "adapt called" << std::endl;
                                    det = Predicates<Precision>::insphere_adapt(permanent, i, sn.vertex(j), simplices,
                                                                                points);
                                }

                                det = Predicates<Precision>::pred_orient(i, simplices, points) *
                                      Predicates<Precision>::pred_insphere(i, sn.vertex(j), simplices, points);

                                if (det < 0) {
//                                    std::cout << "invalid" << std::endl;
                                    ++inv;
                                    valid = false;
                                } else {
                                    ++val;
                                }
                            }
                        }
                    }
                }

                std::cout << "valid: " << val << " invalid: " << inv << std::endl;
            }
        }

        return valid;
    }
};