#pragma once

#include <array>

#include "Datastructures.h"

template<typename Precision>
Precision determinant3x3(const Precision a00, const Precision a01, const Precision a02,
                         const Precision a10, const Precision a11, const Precision a12,
                         const Precision a20, const Precision a21, const Precision a22) {
    return (a00 * (a11 * a22 - a12 * a21) -
            a01 * (a10 * a22 - a12 * a20) +
            a02 * (a10 * a21 - a11 * a20));
}

template<tDimType D, typename Precision, class PointArray>
std::array<Precision, D + 1>
subdeterminants(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                const PointArray &points) {

    Precision ab[D], ac[D], ad[D];
    Precision blift, clift, dlift;

    for (uint i = 0; i < D; ++i) {
        ab[i] = points(pb, i) - points(pa, i);
        ac[i] = points(pc, i) - points(pa, i);
        ad[i] = points(pd, i) - points(pa, i);
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

template<tDimType D, typename Precision, class PointArray>
Precision insphere_fast(const tIndexType &pa,
                        __attribute__((unused)) const tIndexType &pb,
                        __attribute__((unused)) const tIndexType &pc,
                        __attribute__((unused)) const tIndexType &pd,
                        const tIndexType &pe, const std::array<Precision, D + 1> &subdets,
                        const PointArray &points) {

    Precision ae[D];
    Precision elift;

    for (uint i = 0; i < D; ++i) {
        ae[i] = points(pe, i) - points(pa, i);
    }

    elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

    return -ae[X] * subdets[0] + ae[Y] * subdets[1] - ae[Z] * subdets[2] + elift * subdets[3];

}

template<tDimType D, typename Precision, class PointArray>
Precision insphere_fast(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                        const tIndexType &pe, const PointArray &points) {

    Precision ae[D], be[D], ce[D], de[D];
    Precision alift, blift, clift, dlift;
    Precision ab, bc, cd, da, ac, bd;
    Precision abc, bcd, cda, dab;

    for (uint i = 0; i < D; ++i) {
        ae[i] = points(pa, i) - points(pe, i);
        be[i] = points(pb, i) - points(pe, i);
        ce[i] = points(pc, i) - points(pe, i);
        de[i] = points(pd, i) - points(pe, i);
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
}

template<tDimType D, typename Precision, class PointArray, class SimplexArray>
Precision insphere_fast(const tIndexType &pa, const tIndexType &pb, const tIndexType &pc, const tIndexType &pd,
                                                  const tIndexType &pe, const tIndexType &s, const PointArray &points, const
                                         SimplexArray &simplices) {
    
    if constexpr (SimplexArray::hasSubdets){
        return insphere_fast<D, Precision>(pa, pb, pc, pd, pe, simplices.subdets(s), points);
    } else {
        return insphere_fast<D, Precision>(pa, pb, pc, pd, pe, points);
    }
    
}

