#pragma once

#include <array>

#ifdef HAS_Vc

#include <Vc/Vc>

#endif

#include "Datastructures.h"

template<typename Precision>
inline Precision determinant3x3(const Precision a00, const Precision a01, const Precision a02,
                                const Precision a10, const Precision a11, const Precision a12,
                                const Precision a20, const Precision a21, const Precision a22) {
    return (a00 * (a11 * a22 - a12 * a21) -
            a01 * (a10 * a22 - a12 * a20) +
            a02 * (a10 * a21 - a11 * a20));
}

#ifdef HAS_Vc

template<typename Precision>
inline Vc::Vector<Precision>
determinant3x3(const Vc::Vector<Precision> &a00, const Vc::Vector<Precision> &a01, const Vc::Vector<Precision> &a02,
               const Vc::Vector<Precision> &a10, const Vc::Vector<Precision> &a11, const Vc::Vector<Precision> &a12,
               const Vc::Vector<Precision> &a20, const Vc::Vector<Precision> &a21, const Vc::Vector<Precision> &a22) {
    return (a00 * (a11 * a22 - a12 * a21) -
            a01 * (a10 * a22 - a12 * a20) +
            a02 * (a10 * a21 - a11 * a20));
}

#endif

template<tDimType D, typename Precision, class SimplexArray, class PointArray>
void subdeterminants(const tIndexType &s, SimplexArray &simplices, const PointArray &points) {

    if constexpr (SimplexArray::hasSubdets) {
        Precision ab[D], ac[D], ad[D];
        Precision blift, clift, dlift;

        for (uint i = 0; i < D; ++i) {
            ab[i] = points(simplices.vertex(s, 1), i) - points(simplices.vertex(s, 0), i);
            ac[i] = points(simplices.vertex(s, 2), i) - points(simplices.vertex(s, 0), i);
            ad[i] = points(simplices.vertex(s, 3), i) - points(simplices.vertex(s, 0), i);
        }

        blift = ab[X] * ab[X] + ab[Y] * ab[Y] + ab[Z] * ab[Z];
        clift = ac[X] * ac[X] + ac[Y] * ac[Y] + ac[Z] * ac[Z];
        dlift = ad[X] * ad[X] + ad[Y] * ad[Y] + ad[Z] * ad[Z];

        simplices.subdets(s, 0) = determinant3x3<Precision>(ab[Y], ab[Z], blift,
                                                            ac[Y], ac[Z], clift,
                                                            ad[Y], ad[Z], dlift);

        simplices.subdets(s, 1) = determinant3x3<Precision>(ab[X], ab[Z], blift,
                                                            ac[X], ac[Z], clift,
                                                            ad[X], ad[Z], dlift);

        simplices.subdets(s, 2) = determinant3x3<Precision>(ab[X], ab[Y], blift,
                                                            ac[X], ac[Y], clift,
                                                            ad[X], ad[Y], dlift);

        simplices.subdets(s, 3) = determinant3x3<Precision>(ab[X], ab[Y], ab[Z],
                                                            ac[X], ac[Y], ac[Z],
                                                            ad[X], ad[Y], ad[Z]);

    }

}

#ifdef HAS_Vc

template<tDimType D, typename Precision, class SimplexArray, class PointArray>
void subdeterminants(const Vc::Vector<tIndexType> &s, SimplexArray &simplices, const PointArray &points) {

    if constexpr (SimplexArray::hasSubdets) {
        Vc::Vector<Precision> ab[D], ac[D], ad[D];
        Vc::Vector<Precision> blift, clift, dlift;

        for (uint i = 0; i < D; ++i) {
            ab[i] = points(simplices.vertex(s, 1), i) - points(simplices.vertex(s, 0), i);
            ac[i] = points(simplices.vertex(s, 2), i) - points(simplices.vertex(s, 0), i);
            ad[i] = points(simplices.vertex(s, 3), i) - points(simplices.vertex(s, 0), i);
        }

        blift = ab[X] * ab[X] + ab[Y] * ab[Y] + ab[Z] * ab[Z];
        clift = ac[X] * ac[X] + ac[Y] * ac[Y] + ac[Z] * ac[Z];
        dlift = ad[X] * ad[X] + ad[Y] * ad[Y] + ad[Z] * ad[Z];

        simplices.subdets_store(s, 0, determinant3x3<Precision>(ab[Y], ab[Z], blift,
                                                                ac[Y], ac[Z], clift,
                                                                ad[Y], ad[Z], dlift));

        simplices.subdets_store(s, 1, determinant3x3<Precision>(ab[X], ab[Z], blift,
                                                                ac[X], ac[Z], clift,
                                                                ad[X], ad[Z], dlift));

        simplices.subdets_store(s, 2, determinant3x3<Precision>(ab[X], ab[Y], blift,
                                                                ac[X], ac[Y], clift,
                                                                ad[X], ad[Y], dlift));

        simplices.subdets_store(s, 3, determinant3x3<Precision>(ab[X], ab[Y], ab[Z],
                                                                ac[X], ac[Y], ac[Z],
                                                                ad[X], ad[Y], ad[Z]));

    }

}

#endif

template<tDimType D, typename Precision, class SimplexArray, class PointArray>
Precision
insphere_fast(const tIndexType &s, const tIndexType &pe, const SimplexArray &simplices, const PointArray &points) {

    if constexpr (SimplexArray::hasSubdets) {
        Precision ae[D];
        Precision elift;

        for (uint i = 0; i < D; ++i) {
            ae[i] = points(pe, i) - points(simplices.vertex(s, 0), i);
        }

        elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

        return -ae[X] * simplices.subdets(s, 0) + ae[Y] * simplices.subdets(s, 1) - ae[Z] * simplices.subdets(s, 2) +
               elift * simplices.subdets(s, 3);
    } else {
        Precision ae[D], be[D], ce[D], de[D];
        Precision alift, blift, clift, dlift;
        Precision ab, bc, cd, da, ac, bd;
        Precision abc, bcd, cda, dab;

        for (uint i = 0; i < D; ++i) {
            ae[i] = points(simplices.vertex(s, 0), i) - points(pe, i);
            be[i] = points(simplices.vertex(s, 1), i) - points(pe, i);
            ce[i] = points(simplices.vertex(s, 2), i) - points(pe, i);
            de[i] = points(simplices.vertex(s, 3), i) - points(pe, i);
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

}

#ifdef HAS_Vc
template<tDimType D, typename Precision, class SimplexArray, class PointArray>
Vc::Vector<Precision>
insphere_fast(const Vc::Vector<tIndexType> &s, const Vc::Vector<tIndexType> &pe, const SimplexArray &simplices, const PointArray &points) {

    if constexpr (SimplexArray::hasSubdets) {
        Vc::Vector<Precision> ae[D];
        Vc::Vector<Precision> elift;

        for (uint i = 0; i < D; ++i) {
            ae[i] = points(pe, i) - points(simplices.vertex(s, 0), i);
        }

        elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

        return -ae[X] * simplices.subdets(s, 0) + ae[Y] * simplices.subdets(s, 1) - ae[Z] * simplices.subdets(s, 2) +
               elift * simplices.subdets(s, 3);
    } else {
        Vc::Vector<Precision> ae[D], be[D], ce[D], de[D];
        Vc::Vector<Precision> alift, blift, clift, dlift;
        Vc::Vector<Precision> ab, bc, cd, da, ac, bd;
        Vc::Vector<Precision> abc, bcd, cda, dab;

        for (uint i = 0; i < D; ++i) {
            ae[i] = points(simplices.vertex(s, 0), i) - points(pe, i);
            be[i] = points(simplices.vertex(s, 1), i) - points(pe, i);
            ce[i] = points(simplices.vertex(s, 2), i) - points(pe, i);
            de[i] = points(simplices.vertex(s, 3), i) - points(pe, i);
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

}
#endif

