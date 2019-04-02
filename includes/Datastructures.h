#pragma once

#include <vector>
#include <array>

using tIndexType = unsigned long;
constexpr tIndexType INF = ~tIndexType(0);

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
            : m_pointArray(pointArray), m_idx(idx) {

        m_pointArray.ensure(m_idx);
    }

    const Precision &operator[](const tDimType &i) const {
        return m_pointArray(m_idx, i);
    }

    template<typename Ret = Precision &>
    auto operator[](const tDimType &i) -> std::enable_if_t<not std::is_const<PointArray>::value, Ret> {
        return m_pointArray(m_idx, i);
    }
};

template<tDimType D, typename Precision>
class PointAoA {

private:
    std::array<tFloatVector<Precision>, D> coords;

    using tPoint = Point<D, Precision, PointAoA<D, Precision>>;
    using tcPoint = Point<D, Precision, const PointAoA<D, Precision>>;

public:

    Precision insphere_fast(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                            const tIndexType &pe) const {

        Precision aex, bex, cex, dex;
        Precision aey, bey, cey, dey;
        Precision aez, bez, cez, dez;
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        aex = coords[0][pa] - coords[0][pe];
        bex = coords[0][pb] - coords[0][pe];
        cex = coords[0][pc] - coords[0][pe];
        dex = coords[0][pd] - coords[0][pe];
        aey = coords[1][pa] - coords[1][pe];
        bey = coords[1][pb] - coords[1][pe];
        cey = coords[1][pc] - coords[1][pe];
        dey = coords[1][pd] - coords[1][pe];
        aez = coords[2][pa] - coords[2][pe];
        bez = coords[2][pb] - coords[2][pe];
        cez = coords[2][pc] - coords[2][pe];
        dez = coords[2][pd] - coords[2][pe];

        ab = aex * bey - bex * aey;
        bc = bex * cey - cex * bey;
        cd = cex * dey - dex * cey;
        da = dex * aey - aex * dey;

        ac = aex * cey - cex * aey;
        bd = bex * dey - dex * bey;

        abc = aez * bc - bez * ac + cez * ab;
        bcd = bez * cd - cez * bd + dez * bc;
        cda = cez * da + dez * ac + aez * cd;
        dab = dez * ab + aez * bd + bez * da;

        alift = aex * aex + aey * aey + aez * aez;
        blift = bex * bex + bey * bey + bez * bez;
        clift = cex * cex + cey * cey + cez * cez;
        dlift = dex * dex + dey * dey + dez * dez;

        return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
    };

    const Precision &operator()(const tIndexType &i, const tDimType &d) const {
        return coords[d][i];
    }

    Precision &operator()(const tIndexType &i, const tDimType &d) {
        return coords[d][i];
    }

    void ensure(const tIndexType &i) {
        for (tDimType d = 0; d < D; ++d) {
            if (coords[d].size() < i + 1) {
                coords[d].resize(i + 1);
            }
        }
    }

    void ensure(const tIndexType &i) const {
        for (tDimType d = 0; d < D; ++d) {
            if (coords[d].size() < i + 1) {
                throw std::out_of_range(
                        "point index " + std::to_string(i) + " is out of range " + std::to_string(size()));
            }
        }
    }

    tIndexType size() const {
        return coords[0].size();
    }

    tPoint get(const tIndexType &i) {
        return tPoint(*this, i);
    }

    tcPoint get(const tIndexType &i) const {
        return tcPoint(*this, i);
    }

};

template<tDimType D, typename Precision>
class PointPA {

private:
    tFloatVector<Precision> coords;

    using tPoint = Point<D, Precision, PointPA<D, Precision>>;
    using tcPoint = Point<D, Precision, const PointPA<D, Precision>>;

public:

    Precision insphere_fast(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                            const tIndexType &pe) const {

        Precision aex, bex, cex, dex;
        Precision aey, bey, cey, dey;
        Precision aez, bez, cez, dez;
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        aex = coords[pa + 0] - coords[pe + 0];
        bex = coords[pb + 0] - coords[pe + 0];
        cex = coords[pc + 0] - coords[pe + 0];
        dex = coords[pd + 0] - coords[pe + 0];
        aey = coords[pa + 1] - coords[pe + 1];
        bey = coords[pb + 1] - coords[pe + 1];
        cey = coords[pc + 1] - coords[pe + 1];
        dey = coords[pd + 1] - coords[pe + 1];
        aez = coords[pa + 2] - coords[pe + 2];
        bez = coords[pb + 2] - coords[pe + 2];
        cez = coords[pc + 2] - coords[pe + 2];
        dez = coords[pd + 2] - coords[pe + 2];

        ab = aex * bey - bex * aey;
        bc = bex * cey - cex * bey;
        cd = cex * dey - dex * cey;
        da = dex * aey - aex * dey;

        ac = aex * cey - cex * aey;
        bd = bex * dey - dex * bey;

        abc = aez * bc - bez * ac + cez * ab;
        bcd = bez * cd - cez * bd + dez * bc;
        cda = cez * da + dez * ac + aez * cd;
        dab = dez * ab + aez * bd + bez * da;

        alift = aex * aex + aey * aey + aez * aez;
        blift = bex * bex + bey * bey + bez * bez;
        clift = cex * cex + cey * cey + cez * cez;
        dlift = dex * dex + dey * dey + dez * dez;

        return (dlift * abc - clift * dab) + (blift * cda - alift * bcd);
    };

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

    void ensure(const tIndexType &i) const {
        if (coords.size() < D * (i + 1)) {
            throw std::out_of_range(
                    "point index " + std::to_string(i) + " is out of range " + std::to_string(size() / D));
        }
    }

    tIndexType size() const {
        return coords.size() / D;
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
        m_simplexArray.ensure(m_idx);
    }

    const tIndexType &vertex(const tDimType &i) const {
        return m_simplexArray.vertex(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    auto vertex(const tDimType &i) -> std::enable_if_t<not std::is_const<SimplexArray>::value, Ret> {
        return m_simplexArray.vertex(m_idx, i);
    }

    const tIndexType &neighbor(const tDimType &i) const {
        return m_simplexArray.neighbor(m_idx, i);
    }

    template<typename Ret = tIndexType &>
    auto neighbor(const tDimType &i) -> std::enable_if_t<not std::is_const<SimplexArray>::value, Ret> {
        return m_simplexArray.neighbor(m_idx, i);
    }

};

template<tDimType D>
class SimplexAoA {

private:
    std::array<tIndexVector, D + 1> vertices;
    std::array<tIndexVector, D + 1> neighbors;

    using tSimplex = Simplex<D, SimplexAoA<D>>;
    using tcSimplex = Simplex<D, const SimplexAoA<D>>;

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

    void ensure(const tIndexType &i) const {
        for (tDimType d = 0; d < D + 1; ++d) {
            if (std::min(vertices[d].size(), neighbors[d].size()) < i + 1) {
                throw std::out_of_range(
                        "point index " + std::to_string(i) + " is out of range " + std::to_string(size()));
            }
        }
    }

    tIndexType size() const {
        return vertices[0].size();
    }

    tSimplex get(const tIndexType &i) {
        return tSimplex(*this, i);
    }

    tcSimplex get(const tIndexType &i) const {
        return tcSimplex(*this, i);
    }
};

template<tDimType D>
class SimplexPA {

private:
    tIndexVector vertices;
    tIndexVector neighbors;

    using tSimplex = Simplex<D, SimplexPA<D>>;
    using tcSimplex = Simplex<D, const SimplexPA<D>>;

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

    void ensure(const tIndexType &i) const {
        if (std::min(vertices.size(), neighbors.size()) < (D + 1) * (i + 1)) {
            throw std::out_of_range(
                    "point index " + std::to_string(i) + " is out of range " + std::to_string(size() / (D + 1)));
        }
    }

    tIndexType size() const {
        return vertices.size() / (D + 1);
    }

    tSimplex get(const tIndexType &i) {
        return tSimplex(*this, i);
    }

    tcSimplex get(const tIndexType &i) const {
        return tcSimplex(*this, i);
    }
};
