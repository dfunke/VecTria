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

        Precision ae[D], be[D], ce[D], de[D];
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        for (uint d = 0; d < D; ++d) {
            ae[d] = coords[d][pa] - coords[d][pe];
            be[d] = coords[d][pb] - coords[d][pe];
            ce[d] = coords[d][pc] - coords[d][pe];
            de[d] = coords[d][pd] - coords[d][pe];
        }

        ab = ae[0] * be[1] - be[0] * ae[1];
        bc = be[0] * ce[1] - ce[0] * be[1];
        cd = ce[0] * de[1] - de[0] * ce[1];
        da = de[0] * ae[1] - ae[0] * de[1];

        ac = ae[0] * ce[1] - ce[0] * ae[1];
        bd = be[0] * de[1] - de[0] * be[1];

        abc = ae[2] * bc - be[2] * ac + ce[2] * ab;
        bcd = be[2] * cd - ce[2] * bd + de[2] * bc;
        cda = ce[2] * da + de[2] * ac + ae[2] * cd;
        dab = de[2] * ab + ae[2] * bd + be[2] * da;

        alift = ae[0] * ae[0] + ae[1] * ae[1] + ae[2] * ae[2];
        blift = be[0] * be[0] + be[1] * be[1] + be[2] * be[2];
        clift = ce[0] * ce[0] + ce[1] * ce[1] + ce[2] * ce[2];
        dlift = de[0] * de[0] + de[1] * de[1] + de[2] * de[2];

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

        Precision ae[D], be[D], ce[D], de[D];
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        for (uint d = 0; d < D; ++d) {
            ae[d] = coords[pa + d] - coords[pe + d];
            be[d] = coords[pb + d] - coords[pe + d];
            ce[d] = coords[pc + d] - coords[pe + d];
            de[d] = coords[pd + d] - coords[pe + d];
        }

        ab = ae[0] * be[1] - be[0] * ae[1];
        bc = be[0] * ce[1] - ce[0] * be[1];
        cd = ce[0] * de[1] - de[0] * ce[1];
        da = de[0] * ae[1] - ae[0] * de[1];

        ac = ae[0] * ce[1] - ce[0] * ae[1];
        bd = be[0] * de[1] - de[0] * be[1];

        abc = ae[2] * bc - be[2] * ac + ce[2] * ab;
        bcd = be[2] * cd - ce[2] * bd + de[2] * bc;
        cda = ce[2] * da + de[2] * ac + ae[2] * cd;
        dab = de[2] * ab + ae[2] * bd + be[2] * da;

        alift = ae[0] * ae[0] + ae[1] * ae[1] + ae[2] * ae[2];
        blift = be[0] * be[0] + be[1] * be[1] + be[2] * be[2];
        clift = ce[0] * ce[0] + ce[1] * ce[1] + ce[2] * ce[2];
        dlift = de[0] * de[0] + de[1] * de[1] + de[2] * de[2];

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
