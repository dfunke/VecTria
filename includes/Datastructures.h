#pragma once

#include <vector>
#include <array>
#include "MathTools.h"

using tIndexType = unsigned long;
constexpr tIndexType INF = ~tIndexType(0);

using tDimType = unsigned short;

constexpr tDimType X = 0;
constexpr tDimType Y = 1;
constexpr tDimType Z = 2;

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

    std::array<Precision, D + 1>
    subdeterminants(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd) const {

        Precision ab[D], ac[D], ad[D];
        Precision blift, clift, dlift;

        for (uint i = 0; i < D; ++i) {
            ab[i] = coords[i][pb] - coords[i][pa];
            ac[i] = coords[i][pc] - coords[i][pa];
            ad[i] = coords[i][pd] - coords[i][pa];
        }

        blift = ab[X] * ab[X] + ab[Y] * ab[Y] + ab[Z] * ab[Z];
        clift = ac[X] * ac[X] + ac[Y] * ac[Y] + ac[Z] * ac[Z];
        dlift = ad[X] * ad[X] + ad[Y] * ad[Y] + ad[Z] * ad[Z];

        std::array<Precision, D + 1> subdets;

        subdets[0] = determinant3x3<Precision>(ab[Y], ab[Z], blift,
                                               ac[Y], ac[Z], clift,
                                               ad[Y], ad[Z], dlift);

        subdets[1] = determinant3x3<Precision>(ab[X], ab[Z], blift,
                                               ac[X], ac[Z], clift,
                                               ad[X], ad[Z], dlift);

        subdets[2] = determinant3x3<Precision>(ab[X], ab[Y], blift,
                                               ac[X], ac[Y], clift,
                                               ad[X], ad[Y], dlift);

        subdets[3] = determinant3x3<Precision>(ab[X], ab[Y], ab[Z],
                                               ac[X], ac[Y], ac[Z],
                                               ad[X], ad[Y], ad[Z]);

        return subdets;

    }

    Precision insphere_fast(const tIndexType &pa,
                            __attribute__((unused)) const tIndexType &pb,
                            __attribute__((unused)) const tIndexType &pc,
                            __attribute__((unused)) const tIndexType &pd,
                            const tIndexType &pe, const std::array<Precision, D + 1> &subdets) const {

        Precision ae[D];
        Precision elift;

        for (uint i = 0; i < D; ++i) {
            ae[i] = coords[i][pe] - coords[i][pa];
        }

        elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

        return -ae[X] * subdets[0] + ae[Y] * subdets[1] - ae[Z] * subdets[2] + elift * subdets[3];

    }

    Precision insphere_fast(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                            const tIndexType &pe) const {

        Precision ae[D], be[D], ce[D], de[D];
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        for (uint i = 0; i < D; ++i) {
            ae[i] = coords[i][pa] - coords[i][pe];
            be[i] = coords[i][pb] - coords[i][pe];
            ce[i] = coords[i][pc] - coords[i][pe];
            de[i] = coords[i][pd] - coords[i][pe];
        }

        ab = ae[X] * be[Y] - be[X] * ae[Y];
        bc = be[X] * ce[Y] - ce[X] * be[Y];
        cd = ce[X] * de[Y] - de[X] * ce[Y];
        da = de[X] * ae[Y] - ae[X] * de[Y];

        ac = ae[X] * ce[Y] - ce[X] * ae[Y];
        bd = be[X] * de[Y] - de[X] * be[Y];

        abc = ae[Z] * bc - be[Z] * ac + ce[Z] * ab;
        bcd = be[Z] * cd - ce[Z] * bd + de[Z] * bc;
        cda = ce[Z] * da + de[Z] * ac + ae[Z] * cd;
        dab = de[Z] * ab + ae[Z] * bd + be[Z] * da;

        alift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];
        blift = be[X] * be[X] + be[Y] * be[Y] + be[Z] * be[Z];
        clift = ce[X] * ce[X] + ce[Y] * ce[Y] + ce[Z] * ce[Z];
        dlift = de[X] * de[X] + de[Y] * de[Y] + de[Z] * de[Z];

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

    std::array<Precision, D + 1>
    subdeterminants(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd) const {

        Precision ab[D], ac[D], ad[D];
        Precision blift, clift, dlift;

        for (uint i = 0; i < D; ++i) {
            ab[i] = coords[D * pb + i] - coords[D * pa + i];
            ac[i] = coords[D * pc + i] - coords[D * pa + i];
            ad[i] = coords[D * pd + i] - coords[D * pa + i];
        }

        blift = ab[X] * ab[X] + ab[Y] * ab[Y] + ab[Z] * ab[Z];
        clift = ac[X] * ac[X] + ac[Y] * ac[Y] + ac[Z] * ac[Z];
        dlift = ad[X] * ad[X] + ad[Y] * ad[Y] + ad[Z] * ad[Z];

        std::array<Precision, D + 1> subdets;

        subdets[0] = determinant3x3<Precision>(ab[Y], ab[Z], blift,
                                               ac[Y], ac[Z], clift,
                                               ad[Y], ad[Z], dlift);

        subdets[1] = determinant3x3<Precision>(ab[X], ab[Z], blift,
                                               ac[X], ac[Z], clift,
                                               ad[X], ad[Z], dlift);

        subdets[2] = determinant3x3<Precision>(ab[X], ab[Y], blift,
                                               ac[X], ac[Y], clift,
                                               ad[X], ad[Y], dlift);

        subdets[3] = determinant3x3<Precision>(ab[X], ab[Y], ab[Z],
                                               ac[X], ac[Y], ac[Z],
                                               ad[X], ad[Y], ad[Z]);

        return subdets;

    }

    Precision insphere_fast(const tIndexType &pa,
                            __attribute__((unused)) const tIndexType &pb,
                            __attribute__((unused)) const tIndexType &pc,
                            __attribute__((unused)) const tIndexType &pd,
                            const tIndexType &pe, const std::array<Precision, D + 1> &subdets) const {

        Precision ae[D];
        Precision elift;

        for (uint i = 0; i < D; ++i) {
            ae[i] = coords[D * pe + i] - coords[D * pa + i];
        }

        elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

        return -ae[X] * subdets[0] + ae[Y] * subdets[1] - ae[Z] * subdets[2] + elift * subdets[3];

    }

    Precision insphere_fast(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                            const tIndexType &pe) const {

        Precision ae[D], be[D], ce[D], de[D];
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        for (uint i = 0; i < D; ++i) {
            ae[i] = coords[D * pa + i] - coords[D * pe + i];
            be[i] = coords[D * pb + i] - coords[D * pe + i];
            ce[i] = coords[D * pc + i] - coords[D * pe + i];
            de[i] = coords[D * pd + i] - coords[D * pe + i];
        }

        ab = ae[X] * be[Y] - be[X] * ae[Y];
        bc = be[X] * ce[Y] - ce[X] * be[Y];
        cd = ce[X] * de[Y] - de[X] * ce[Y];
        da = de[X] * ae[Y] - ae[X] * de[Y];

        ac = ae[X] * ce[Y] - ce[X] * ae[Y];
        bd = be[X] * de[Y] - de[X] * be[Y];

        abc = ae[Z] * bc - be[Z] * ac + ce[Z] * ab;
        bcd = be[Z] * cd - ce[Z] * bd + de[Z] * bc;
        cda = ce[Z] * da + de[Z] * ac + ae[Z] * cd;
        dab = de[Z] * ab + ae[Z] * bd + be[Z] * da;

        alift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];
        blift = be[X] * be[X] + be[Y] * be[Y] + be[Z] * be[Z];
        clift = ce[X] * ce[X] + ce[Y] * ce[Y] + ce[Z] * ce[Z];
        dlift = de[X] * de[X] + de[Y] * de[Y] + de[Z] * de[Z];

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
