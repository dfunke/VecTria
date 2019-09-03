#pragma once

#include <vector>
#include <array>

#include "Datastructures.h"
#include "MathTools.h"

template<tDimType pD,
        typename pPrecision,
        template<typename, tDimType> class pMemoryLayout,
        template<class> class pPrecomputeStrategy,
        template<class> class pOppVertex>
struct Traits {

    template<typename T, tDimType D2>
    using MemoryLayout = pMemoryLayout<T, D2>;

    template<class T>
    using PrecomputeStrategy = pPrecomputeStrategy<Traits>;

    template<class T>
    using OppVertex = pOppVertex<Traits>;

    using Precision = pPrecision;

    static constexpr tDimType D = pD;
};

template<tDimType D, typename Precision, class PointArray>
class Point {

private:

    PointArray &m_pointArray;
    tIndexType m_idx;

public:

    Point(PointArray &pointArray, const tIndexType &idx)
            : m_pointArray(pointArray), m_idx(idx) {

        if constexpr (not std::is_const_v<PointArray>) {
            m_pointArray.ensure(m_idx);
        }
    }

    inline const Precision &operator[](const tDimType &i) const {
        return m_pointArray.coords(m_idx, i);
    }

    template<typename Ret = Precision &>
    inline auto operator[](const tDimType &i) -> std::enable_if_t<not std::is_const_v<PointArray>, Ret> {
        return m_pointArray.coords(m_idx, i);
    }

    template<class OtherPoint>
    bool operator==(const OtherPoint &p) const {

        for (tDimType d = 0; d < D; ++d) {
            if (operator[](d) != p[d])
                return false;
        }

        return true;
    }
};

template<class Traits>
class PointArray {

public:
    template<typename T, tDimType D2>
    using MemoryLayout = typename Traits::template MemoryLayout<T, D2>;

    using Precision = typename Traits::Precision;
    static constexpr tDimType D = Traits::D;
    static constexpr bool isVectorized = MemoryLayout<tIndexType, D + 1>::isVectorized;

public:
    MemoryLayout<Precision, D> coords;

private:
    using tPoint = Point<D, Precision, PointArray>;
    using tcPoint = Point<D, Precision, const PointArray>;

public:

    void ensure(const tIndexType &i) {
        coords.ensure(i);
    }

    void ensure(const tIndexType &i) const {
        coords.ensure(i);
    }

    tIndexType size() const {
        return coords.size();
    }

    tPoint get(const tIndexType &i) {
        return tPoint(*this, i);
    }

    tcPoint get(const tIndexType &i) const {
        return tcPoint(*this, i);
    }

    MemoryLayout<Precision, D> getDelta(const PointArray &b) const {
        assert(coords.size() == b.size());

        MemoryLayout<Precision, D> delta;
        delta.ensure(coords.size());

        if constexpr (isVectorized) {
            for (tIndexType i = 0;
                 i + Vc::Vector<tIndexType>::size() - 1 < coords.size(); i += Vc::Vector<tIndexType>::size()) {
                for (tDimType d = 0; d < D; ++d) {
                    delta.store(i, d, b.coords.vec(i, d) - coords.vec(i, d));
                }
            }
        } else {
            for (tIndexType i = 0; i < coords.size(); ++i) {
                for (tDimType d = 0; d < D; ++d) {
                    delta(i, d) = b.coords(i, d) - coords(i, d);
                }
            }
        }

        return delta;
    }

};

template<tDimType D, class SimplexArray>
class Simplex {

private:

    SimplexArray &m_simplexArray;
    tIndexType m_idx;

public:

    Simplex(SimplexArray &simplexArray, const tIndexType &idx)
            : m_simplexArray(simplexArray), m_idx(idx) {

        if constexpr (not std::is_const_v<SimplexArray>) {
            m_simplexArray.ensure(m_idx);
        }
    }

    inline const tIndexType &vertex(const tDimType &i) const {
        return m_simplexArray.vertices(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    inline auto vertex(const tDimType &i) -> std::enable_if_t<not std::is_const_v<SimplexArray>, Ret> {
        return m_simplexArray.vertices(m_idx, i);
    }

    inline const tIndexType &neighbor(const tDimType &i) const {
        return m_simplexArray.neighbors(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    inline auto neighbor(const tDimType &i) -> std::enable_if_t<not std::is_const_v<SimplexArray>, Ret> {
        return m_simplexArray.neighbors(m_idx, i);
    }

    template<typename Ret = const tIndexType &>
    inline auto oppVertex(const tDimType &i) const -> std::enable_if_t<SimplexArray::hasOppVertex, Ret> {
        return m_simplexArray.opp_vertex(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    inline auto oppVertex(const tDimType &i) -> std::enable_if_t<
            SimplexArray::hasOppVertex and not std::is_const_v<SimplexArray>, Ret> {
        return m_simplexArray.opp_vertex(m_idx, i);
    }

    template<class OtherSimplex>
    bool operator==(const OtherSimplex &s) const {

        for (tDimType d = 0; d < D + 1; ++d) {
            if (vertex(d) != s.vertex(d) || neighbor(d) != s.neighbor(d))
                return false;
        }

        return true;
    }

};

template<class Traits>
class NoPrecomputation {

public:

    static constexpr bool hasSubdets = false;

    template<class SimplexArray, class PointArray>
    void precompute(__attribute__((unused)) SimplexArray &simplices,
                    __attribute__((unused)) const PointArray &points) {}
};

template<class Traits>
class PrecomputeSubDets {

public:
    template<typename T, tDimType D2>
    using MemoryLayout = typename Traits::template MemoryLayout<T, D2>;

    using Precision = typename Traits::Precision;
    static constexpr tDimType D = Traits::D;

public:
    MemoryLayout<Precision, D + 1> subdets;

public:

    static constexpr bool hasSubdets = true;

    template<class SimplexArray, class PointArray>
    void precompute(SimplexArray &simplices, const PointArray &points) {

        subdets.ensure(simplices.size() - 1);

        if constexpr (MemoryLayout<Precision, D + 1>::isVectorized) {
            for (tIndexType i = 0;
                 i + Vc::Vector<tIndexType>::size() - 1 < simplices.size(); i += Vc::Vector<tIndexType>::size()) {

                Predicates<Precision>::subdeterminants_v(i, simplices, points);

            }
        } else {

            for (tIndexType i = 0; i < simplices.size(); ++i) {
                Predicates<Precision>::subdeterminants(i, simplices, points);

            }
        }

    }
};

template<class Traits>
class NoOppVertex {

public:

    static constexpr bool hasOppVertex = false;

};

template<class Traits>
class WithOppVertex {

public:
    template<typename T, tDimType D2>
    using MemoryLayout = typename Traits::template MemoryLayout<T, D2>;

    static constexpr tDimType D = Traits::D;
    static constexpr bool hasOppVertex = true;

public:
    MemoryLayout<tIndexType, D + 1> opp_vertex;
};

template<class Traits>
class SimplexArray : public Traits::template PrecomputeStrategy<Traits>, public Traits::template OppVertex<Traits> {

public:
    template<typename T, tDimType D2>
    using MemoryLayout = typename Traits::template MemoryLayout<T, D2>;

    using Precision = typename Traits::Precision;
    static constexpr tDimType D = Traits::D;
    static constexpr bool isVectorized = MemoryLayout<tIndexType, D + 1>::isVectorized;
    static constexpr bool hasOppVertex = Traits::template OppVertex<Traits>::hasOppVertex;

public:
    MemoryLayout<tIndexType, D + 1> vertices;
    MemoryLayout<tIndexType, D + 1> neighbors;

private:

    using pcBase = typename Traits::template PrecomputeStrategy<Traits>;
    using ovBase = typename Traits::template OppVertex<Traits>;
    using tSimplex = Simplex<D, SimplexArray>;
    using tcSimplex = Simplex<D, const SimplexArray>;

public:

    void ensure(const tIndexType &i) {
        vertices.ensure(i);
        neighbors.ensure(i);

        if constexpr (hasOppVertex) {
            ovBase::opp_vertex.ensure(i);
        }
    }

    void ensure(const tIndexType &i) const {
        vertices.ensure(i);
        neighbors.ensure(i);

        if constexpr (hasOppVertex) {
            ovBase::opp_vertex.ensure(i);
        }
    }

    tIndexType size() const {
        return vertices.size();
    }

    tSimplex get(const tIndexType &i) {
        return tSimplex(*this, i);
    }

    tcSimplex get(const tIndexType &i) const {
        return tcSimplex(*this, i);
    }

    template<class PointArray>
    void precompute(const PointArray &points) {
        pcBase::precompute(*this, points);
    }
};
