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

template<tDimType D, typename Precision>
using tVector = std::array<Precision, D>;

template<typename Vector>
struct Box {
    
    Box() { }
    
    Box(const Vector &_low, const Vector &_high) : low(_low), high(_high) { }
    
    Vector low;
    Vector high;
    
    template<typename Point>
    bool contains(const Point & p) const {
        static_assert(Point::D == std::tuple_size<Vector>::value);
        
        for(tDimType d = 0; d < Point::D; ++d){
            if(!(low[d] <= p[d] && p[d] < high[d]))
                return false;
        }
        
        return true;
    }
};

template<class PointArray>
class Point {

private:

    PointArray &m_pointArray;
    tIndexType m_idx;

public:
    using Precision = typename PointArray::Precision;
    static constexpr tDimType D = PointArray::D;

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
    using tPoint = Point<PointArray>;
    using tcPoint = Point<const PointArray>;

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

template<class SimplexArray>
class Simplex {

private:

    SimplexArray &m_simplexArray;
    tIndexType m_idx;

public:
    static constexpr tDimType D = SimplexArray::D;

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
    using tSimplex = Simplex<SimplexArray>;
    using tcSimplex = Simplex<const SimplexArray>;

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

template<class PointArray, class OStream>
OStream &operator<<(OStream &os, const Point<PointArray> &p) {
    os << "[" << p[0];

    for (tDimType d = 1; d < PointArray::D; ++d) {
        os << ", " << p[d];

    }
    os << "]";

    return os;
}

template<class SimplexArray, class OStream>
OStream &operator<<(OStream &os, const Simplex<SimplexArray> &s) {
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

template<class Point, class OStream>
OStream &operator<<(OStream &os, const Box<Point> &b) {
    os << "[" << b.low[0];

    for (tDimType d = 1; d < Box<Point>::D; ++d) {
        os << ", " << b.low[d];

    }
    os << "]";
    os << " [" << b.high[0];

    for (tDimType d = 1; d < Box<Point>::D; ++d) {
        os << ", " << b.high[d];

    }
    os << "]";

    return os;
}
