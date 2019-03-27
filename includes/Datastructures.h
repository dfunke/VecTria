#pragma once

#include <vector>
#include <array>

using tIndexType = unsigned long;
using tDimType = unsigned short;

template<typename Precision>
using tFloatVector = std::vector<Precision>;

using tIndexVector = std::vector<tIndexType>;

template<tDimType D, typename Precision, class PointArray>
class Point {

private:

    PointArray &m_pointArray;
    tIndexType m_idx;

public:

    Point(PointArray &pointArray, const tIndexType &idx)
            : m_pointArray(pointArray), m_idx(idx) {}

    const Precision &operator[](const tDimType &i) const {
        return m_pointArray(m_idx, i);
    }

    Precision &operator[](const tDimType &i) {
        return m_pointArray(m_idx, i);
    }
};

template<tDimType D, typename Precision>
class PointAoA {

private:
    std::array<tFloatVector<Precision>, D> coords;

    using tPoint = Point<D, Precision, PointAoA<D, Precision>>;

public:

    const Precision &operator()(const tIndexType &i, const tDimType &d) const {
        return coords[d][i];
    }

    Precision &operator()(const tIndexType &i, const tDimType &d) {
        return coords[d][i];
    }

    void ensure(const tIndexType &i) {
        for (tDimType d = 0; d < D + 1; ++d) {
            if (coords[d].size() < i + 1) {
                coords[d].resize(i + 1);
            }
        }
    }

    tPoint get(const tIndexType &i) {
        return tPoint(*this, i);
    }

};

template<tDimType D, typename Precision>
class PointPA {

private:
    tFloatVector<Precision> coords;

    using tPoint = Point<D, Precision, PointPA<D, Precision>>;

public:

    const Precision &operator()(const tIndexType &i, const tDimType &d) const {
        return coords[D * i + d];
    }

    Precision &operator()(const tIndexType &i, const tDimType &d) {
        return coords[D * i + d];
    }

    void ensure(const tIndexType &i) {
        if (coords.size() < D * (i + 1)) {
            coords.resize(D * (i + 1));
        }
    }

    tPoint get(const tIndexType &i) {
        return tPoint(*this, i);
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
        m_simplexArray.ensure(m_idx);
    }

    const tIndexType &vertex(const tDimType &i) const {
        return m_simplexArray.vertices[i][m_idx];
    }

    tIndexType &vertex(const tDimType &i) {
        return m_simplexArray.vertices[i][m_idx];
    }

    const tIndexType &neighbor(const tDimType &i) const {
        return m_simplexArray.neighbors[i][m_idx];
    }

    tIndexType &neighbor(const tDimType &i) {
        return m_simplexArray.neighbors[i][m_idx];
    }

};

template<tDimType D>
class SimplexAoA {

private:
    std::array<tIndexVector, D + 1> vertices;
    std::array<tIndexVector, D + 1> neighbors;

    using tSimplex = Simplex<D, SimplexAoA<D>>;

public:

    const tIndexType &vertex(const tIndexType &i, const tDimType &d) const {
        return vertices[d][i];
    }

    tIndexType &vertex(const tIndexType &i, const tDimType &d) {
        return vertices[d][i];
    }

    const tIndexType &neighbor(const tIndexType &i, const tDimType &d) const {
        return neighbors[d][i];
    }

    tIndexType &neighbor(const tIndexType &i, const tDimType &d) {
        return neighbors[d][i];
    }

    void ensure(const tIndexType &i) {
        for (tDimType d = 0; d < D + 1; ++d) {
            if (vertices[d].size() < i + 1) {
                vertices[d].resize(i + 1);
            }

            if (neighbors[d].size() < i + 1) {
                neighbors[d].resize(i + 1);
            }
        }
    }

    tSimplex get(const tIndexType &i) {
        return tSimplex(*this, i);
    }
};

template<tDimType D>
class SimplexPA {

private:
    tIndexVector vertices;
    tIndexVector neighbors;

    using tSimplex = Simplex<D, SimplexPA<D>>;

public:

    const tIndexType &vertex(const tIndexType &i, const tDimType &d) const {
        return vertices[(D + 1) * i + d];
    }

    tIndexType &vertex(const tIndexType &i, const tDimType &d) {
        return vertices[(D + 1) * i + d];
    }

    const tIndexType &neighbor(const tIndexType &i, const tDimType &d) const {
        return neighbors[(D + 1) * i + d];
    }

    tIndexType &neighbor(const tIndexType &i, const tDimType &d) {
        return neighbors[(D + 1) * i + d];
    }

    void ensure(const tIndexType &i) {
        if (vertices.size() < (D + 1) * (i + 1)) {
            vertices.resize((D + 1) * (i + 1));
        }

        if (neighbors.size() < (D + 1) * (i + 1)) {
            neighbors.resize((D + 1) * (i + 1));
        }
    }

    tSimplex get(const tIndexType &i) {
        return tSimplex(*this, i);
    }
};
