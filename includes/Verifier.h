#pragma once

#include <ostream>

#include "GeometryStructures.h"

template<class PointArray>
std::ostream &operator<<(std::ostream &os, const Point<PointArray> &p) {
    os << "[" << p[0];

    for (tDimType d = 1; d < PointArray::D; ++d) {
        os << ", " << p[d];

    }
    os << "]";

    return os;
}

template<class SimplexArray>
std::ostream &operator<<(std::ostream &os, const Simplex<SimplexArray> &s) {
    os << "[" << s.vertex(0);

    for (tDimType d = 1; d < SimplexArray::D + 1; ++d) {
        os << ", " << s.vertex(d);

    }
    os << "]";
    os << " (" << s.neighbor(0);

    for (tDimType d = 1; d < SimplexArray::D + 1; ++d) {
        os << ", " << s.neighbor(d);

    }
    os << ")";

    return os;
}

template<class OutStream>
struct Verifier {
    Verifier(OutStream& outStream) : outStream(outStream) {
    }

    template<class SimplexArray1, class PointArray1, class SimplexArray2, class PointArray2>
    bool verify(const SimplexArray1 &simplices1, const PointArray1 &points1,
                const SimplexArray2 &simplices2, const PointArray2 &points2) {
        outStream << "Checking\n"
                  << "\tLayout1: "
                  << SimplexArray1::template MemoryLayout<typename SimplexArray1::Precision, SimplexArray1::D>::name()
                  << "\tLayout2: "
                  << SimplexArray2::template MemoryLayout<typename SimplexArray2::Precision, SimplexArray2::D>::name()
                  << std::endl;

        // check number of points
        if (points1.size() != points2.size()) {
            outStream << "Number points: " << points1.size() << " vs " << points2.size() << std::endl;
            return false;
        }

        // check points
        for (tIndexType i = 0; i < points1.size(); ++i) {
            auto p1 = points1.get(i);
            auto p2 = points2.get(i);

            if (!(p1 == p2)) {
                outStream << "Points " << i << " differ: (" << p1 << ") vs (" << p2 << ")" << std::endl;
                return false;
            }
        }

        // check number of simplices
        if (simplices1.size() != simplices2.size()) {
            outStream << "Number simplices: " << simplices1.size() << " vs " << simplices2.size() << std::endl;
            return false;
        }

        // check simplices
        for (tIndexType i = 0; i < simplices1.size(); ++i) {
            auto s1 = simplices1.get(i);
            auto s2 = simplices2.get(i);

            if (!(s1 == s2)) {
                outStream << "Simplices " << i << " differ: (" << s1 << ") vs (" << s2 << ")" << std::endl;
                return false;
            }
        }

        bool result = true;
        if constexpr (SimplexArray2::isVectorized) {
            // vector of indices access
            for (Vc::Vector<tIndexType> i = Vc::Vector<tIndexType>::IndexesFromZero();
                 i.max() < simplices2.size(); i += Vc::Vector<tIndexType>(Vc::Vector<tIndexType>::size())) {


                Vc::Vector<tIndexType> v[SimplexArray2::D + 1], n[SimplexArray2::D + 1];

                for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                    v[d] = simplices2.vertices.vec(i, d);
                    n[d] = simplices2.neighbors.vec(i, d);
                }

                for (std::size_t j = 0; j < Vc::Vector<tIndexType>::size(); ++j) {
                    auto s1 = simplices1.get(i[j]);
                    bool valid = true;

                    for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                        if ((s1.vertex(d) != v[d][j]) || (s1.neighbor(d) != n[d][j])) {
                            valid = false;
                        }
                    }

                    if (!valid) {
                        outStream << "Simplices " << j << ": " << i[j]
                                  << " differ: (" << s1 << ") vs (" << "[" << v[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            outStream << ", " << v[d][j];

                        }
                        outStream << "]";
                        outStream << " (" << n[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            outStream << ", " << n[d][j];

                        }
                        outStream << ")" << ")" << std::endl;
                        result = false;
                    }
                }
            }

            // contiguous load access
            for (tIndexType i = 0;
                 i + Vc::Vector<tIndexType>::size() - 1 < simplices2.size(); i += Vc::Vector<tIndexType>::size()) {

                Vc::Vector<tIndexType> v[SimplexArray2::D + 1], n[SimplexArray2::D + 1];

                for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                    v[d] = simplices2.vertices.vec(i, d);
                    n[d] = simplices2.neighbors.vec(i, d);
                }

                for (std::size_t j = 0; j < Vc::Vector<tIndexType>::size(); ++j) {
                    auto s1 = simplices1.get(i + j);
                    bool valid = true;

                    for (tDimType d = 0; d < SimplexArray2::D + 1; ++d) {
                        if ((s1.vertex(d) != v[d][j]) || (s1.neighbor(d) != n[d][j])) {
                            valid = false;
                        }
                    }

                    if (!valid) {
                        outStream << "Simplices " << j << ": " << i + j
                                  << " differ: (" << s1 << ") vs (" << "[" << v[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            outStream << ", " << v[d][j];

                        }
                        outStream << "]";
                        outStream << " (" << n[0][j];

                        for (tDimType d = 1; d < SimplexArray2::D + 1; ++d) {
                            outStream << ", " << n[d][j];

                        }
                        outStream << ")" << ")" << std::endl;
                        result = false;
                    }
                }
            }
        }
        return result;
    }

private:
    OutStream& outStream;
};