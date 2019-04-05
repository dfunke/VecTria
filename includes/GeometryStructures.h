#pragma once

#include <vector>
#include <array>

#include "Datastructures.h"
#include "MathTools.h"

template<tDimType D, typename Precision, class PointArray>
class Point {

private:

    PointArray &m_pointArray;
    tIndexType m_idx;

public:

    Point(PointArray &pointArray, const tIndexType &idx)
            : m_pointArray(pointArray), m_idx(idx) {

        if constexpr (not std::is_const_v<PointArray>){
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

template<tDimType D, typename Precision, template<typename, tDimType> class MemoryLayout>
class PointArray {

private:
    MemoryLayout<Precision, D> coords;

    using tPoint = Point<D, Precision, PointArray<D, Precision, MemoryLayout>>;
    using tcPoint = Point<D, Precision, const PointArray<D, Precision, MemoryLayout>>;

public:

    inline const Precision &operator()(const tIndexType &i, const tDimType &d) const {
        return coords(i, d);
    }

    inline Precision &operator()(const tIndexType &i, const tDimType &d) {
        return coords(i, d);
    }

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
        
        if constexpr (not std::is_const_v<SimplexArray>){
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

template<class SimplexArray>
class NoPrecomputation {
    
public:
    template<class PointArray>
    void precompute(const SimplexArray &simplices, const PointArray &points){ }
};

template<class SimplexArray>
class PrecomputeSubDets {

private:
    typename SimplexArray::template tMemoryLayout<SimplexArray::tPrecision, SimplexArray::cD + 1> subdets;
    
public:
    
    template<class PointArray>
    void precompute(const SimplexArray &simplices, const PointArray &points){
        for(tIndexType i = 0; i < simplices.size(); ++i){
            auto s = simplices.get(i);
            auto subdet = subdeterminants(s.vertices(0), s.vertices(1), s.vertices(2), s.vertices(3), points);
            
            for(tDimType d = 0; d < SimplexArray::cD + 1; ++d){
                subdets(i, d) = subdet[d];
            }
        }
    }
};

template<tDimType D, typename Precision, template<typename, tDimType> class MemoryLayout, template<class> class PrecomputeStrategy>
class SimplexArray {

public:
    template<typename T, tDimType D2>
    using tMemoryLayout = MemoryLayout<T, D2>;
    
    using tPrecision = Precision;
    static constexpr tDimType cD = D;
    
private:
    MemoryLayout<tIndexType, D + 1> vertices;
    MemoryLayout<tIndexType, D + 1> neighbors;
    
    PrecomputeStrategy<SimplexArray> precomputeStrategy;

    using tSimplex = Simplex<D, SimplexArray>;
    using tcSimplex = Simplex<D, const SimplexArray>;

public:

    inline const tIndexType &vertex(const tIndexType &i, const tDimType &d) const {
        return vertices(i, d);
    }

    inline tIndexType &vertex(const tIndexType &i, const tDimType &d) {
        return vertices(i, d);
    }

    inline const tIndexType &neighbor(const tIndexType &i, const tDimType &d) const {
        return neighbors(i, d);
    }

    inline tIndexType &neighbor(const tIndexType &i, const tDimType &d) {
        return neighbors(i, d);
    }

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
        precomputeStrategy.precompute(*this, points);
    }
};
