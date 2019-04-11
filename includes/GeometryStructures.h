#pragma once

#include <vector>
#include <array>

#include "Datastructures.h"
#include "MathTools.h"

template<tDimType pD,
        typename pPrecision,
        template<typename, tDimType> class pMemoryLayout,
        template<class> class pPrecomputeStrategy>
struct Traits {

    template<typename T, tDimType D2>
    using MemoryLayout = pMemoryLayout<T, D2>;

    template<class T>
    using PrecomputeStrategy = pPrecomputeStrategy<Traits>;

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
        return m_pointArray(m_idx, i);
    }

    template<typename Ret = Precision &>
    inline auto operator[](const tDimType &i) -> std::enable_if_t<not std::is_const_v<PointArray>, Ret> {
        return m_pointArray(m_idx, i);
    }
};

template<class Traits>
class PointArray {

public:
    template<typename T, tDimType D2>
    using MemoryLayout = typename Traits::template MemoryLayout<T, D2>;

    using Precision = typename Traits::Precision;
    static constexpr tDimType D = Traits::D;

private:
    MemoryLayout<Precision, D> coords;

    using tPoint = Point<D, Precision, PointArray>;
    using tcPoint = Point<D, Precision, const PointArray>;

public:

    inline const Precision &operator()(const tIndexType &i, const tDimType &d) const {
        return coords(i, d);
    }

    inline Precision &operator()(const tIndexType &i, const tDimType &d) {
        return coords(i, d);
    }

#ifdef HAS_Vc

    inline Vc::Vector<Precision> operator()(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        return coords(i, d);
    }

#endif

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
        return m_simplexArray.vertex(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    inline auto vertex(const tDimType &i) -> std::enable_if_t<not std::is_const_v<SimplexArray>, Ret> {
        return m_simplexArray.vertex(m_idx, i);
    }

    inline const tIndexType &neighbor(const tDimType &i) const {
        return m_simplexArray.neighbor(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    inline auto neighbor(const tDimType &i) -> std::enable_if_t<not std::is_const_v<SimplexArray>, Ret> {
        return m_simplexArray.neighbor(m_idx, i);
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

private:
    MemoryLayout<Precision, D + 1> f_subdets;

public:

    static constexpr bool hasSubdets = true;

    template<class SimplexArray, class PointArray>
    void precompute(SimplexArray &simplices, const PointArray &points) {

        f_subdets.ensure(simplices.size() - 1);

#ifdef HAS_Vc
        for (auto i = Vc::Vector<tIndexType>::IndexesFromZero(); i.max() < simplices.size();
             i += Vc::Vector<tIndexType>(Vc::Vector<tIndexType>::size())) {

            subdeterminants<D, Precision>(i, simplices, points);

        }
#else
        for (tIndexType i = 0; i < simplices.size(); ++i) {
            subdeterminants<D, Precision>(i, simplices, points);

        }
#endif
    }

    inline const Precision &subdets(const tIndexType &i, const tDimType &d) const {
        return f_subdets(i, d);
    }

    inline Precision &subdets(const tIndexType &i, const tDimType &d) {
        return f_subdets(i, d);
    }

#ifdef HAS_Vc

    inline Vc::Vector<Precision> subdets(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        return f_subdets(i, d);
    }

#endif
};

template<class Traits>
class SimplexArray : public Traits::template PrecomputeStrategy<Traits> {

public:
    template<typename T, tDimType D2>
    using MemoryLayout = typename Traits::template MemoryLayout<T, D2>;

    using Precision = typename Traits::Precision;
    static constexpr tDimType D = Traits::D;

private:
    MemoryLayout<tIndexType, D + 1> vertices;
    MemoryLayout<tIndexType, D + 1> neighbors;

    using base = typename Traits::template PrecomputeStrategy<Traits>;
    using tSimplex = Simplex<D, SimplexArray>;
    using tcSimplex = Simplex<D, const SimplexArray>;

public:

    inline const tIndexType &vertex(const tIndexType &i, const tDimType &d) const {
        return vertices(i, d);
    }

    inline tIndexType &vertex(const tIndexType &i, const tDimType &d) {
        return vertices(i, d);
    }

#ifdef HAS_Vc

    inline Vc::Vector<tIndexType> vertex(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        return vertices(i, d);
    }

#endif

    inline const tIndexType &neighbor(const tIndexType &i, const tDimType &d) const {
        return neighbors(i, d);
    }

    inline tIndexType &neighbor(const tIndexType &i, const tDimType &d) {
        return neighbors(i, d);
    }

#ifdef HAS_Vc

    inline Vc::Vector<tIndexType> neighbor(const Vc::Vector<tIndexType> &i, const tDimType &d) const {
        return neighbor(i, d);
    }

#endif

    void ensure(const tIndexType &i) {
        vertices.ensure(i);
        neighbors.ensure(i);
    }

    void ensure(const tIndexType &i) const {
        vertices.ensure(i);
        neighbors.ensure(i);
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
        base::precompute(*this, points);
    }
};
