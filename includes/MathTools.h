#pragma once

#include <array>

#ifdef HAS_Vc

#include <Vc/Vc>

#endif

#include "Datastructures.h"
#include "Predicates.h"

template<typename Precision>
class Predicates : public PredicatesBase<Precision> {

private:
    using super = PredicatesBase<Precision>;

public:

    static inline Precision determinant3x3(const Precision a00, const Precision a01, const Precision a02,
                                           const Precision a10, const Precision a11, const Precision a12,
                                           const Precision a20, const Precision a21, const Precision a22) {
        return (a00 * (a11 * a22 - a12 * a21) -
                a01 * (a10 * a22 - a12 * a20) +
                a02 * (a10 * a21 - a11 * a20));
    }

#ifdef HAS_Vc

    static inline Vc::Vector<Precision>
    determinant3x3(const Vc::Vector<Precision> &a00, const Vc::Vector<Precision> &a01, const Vc::Vector<Precision> &a02,
                   const Vc::Vector<Precision> &a10, const Vc::Vector<Precision> &a11, const Vc::Vector<Precision> &a12,
                   const Vc::Vector<Precision> &a20, const Vc::Vector<Precision> &a21,
                   const Vc::Vector<Precision> &a22) {
        return (a00 * (a11 * a22 - a12 * a21) -
                a01 * (a10 * a22 - a12 * a20) +
                a02 * (a10 * a21 - a11 * a20));
    }

#endif

    static bool isFinal(const Precision det, const Precision permanent) {
        Precision errbound = super::isperrboundA * permanent;
        return ((det > errbound) || (-det > errbound));
    }

#ifdef HAS_Vc

    static inline auto isFinal(const Vc::Vector<Precision> &det, const Vc::Vector<Precision> &permanent) {
        Vc::Vector<Precision> errbound = super::isperrboundA * permanent;
        return ((det > errbound) || (-det > errbound));
    }

#endif

    template<class SimplexArray, class PointArray>
    static void subdeterminants(const tIndexType &s, SimplexArray &simplices, const PointArray &points) {

        const tDimType D = SimplexArray::D;

        if constexpr (SimplexArray::hasSubdets) {
            Precision ab[D], ac[D], ad[D];
            Precision blift, clift, dlift;

            for (uint i = 0; i < D; ++i) {
                ab[i] = points.coords(simplices.vertices(s, 1), i) - points.coords(simplices.vertices(s, 0), i);
                ac[i] = points.coords(simplices.vertices(s, 2), i) - points.coords(simplices.vertices(s, 0), i);
                ad[i] = points.coords(simplices.vertices(s, 3), i) - points.coords(simplices.vertices(s, 0), i);
            }

            blift = ab[X] * ab[X] + ab[Y] * ab[Y] + ab[Z] * ab[Z];
            clift = ac[X] * ac[X] + ac[Y] * ac[Y] + ac[Z] * ac[Z];
            dlift = ad[X] * ad[X] + ad[Y] * ad[Y] + ad[Z] * ad[Z];

            simplices.subdets(s, 0) = determinant3x3(ab[Y], ab[Z], blift,
                                                     ac[Y], ac[Z], clift,
                                                     ad[Y], ad[Z], dlift);

            simplices.subdets(s, 1) = determinant3x3(ab[X], ab[Z], blift,
                                                     ac[X], ac[Z], clift,
                                                     ad[X], ad[Z], dlift);

            simplices.subdets(s, 2) = determinant3x3(ab[X], ab[Y], blift,
                                                     ac[X], ac[Y], clift,
                                                     ad[X], ad[Y], dlift);

            simplices.subdets(s, 3) = determinant3x3(ab[X], ab[Y], ab[Z],
                                                     ac[X], ac[Y], ac[Z],
                                                     ad[X], ad[Y], ad[Z]);

        }

    }

#ifdef HAS_Vc

    template<class SimplexArray, class PointArray, typename Simplex>
    static void subdeterminants_v(const Simplex &s, SimplexArray &simplices, const PointArray &points) {

        const tDimType D = SimplexArray::D;

        if constexpr (SimplexArray::hasSubdets) {
            Vc::Vector<Precision> ab[D], ac[D], ad[D];
            Vc::Vector<Precision> blift, clift, dlift;

            for (uint i = 0; i < D; ++i) {
                ab[i] = points.coords.vec(simplices.vertices.vec(s, 1), i) -
                        points.coords.vec(simplices.vertices.vec(s, 0), i);
                ac[i] = points.coords.vec(simplices.vertices.vec(s, 2), i) -
                        points.coords.vec(simplices.vertices.vec(s, 0), i);
                ad[i] = points.coords.vec(simplices.vertices.vec(s, 3), i) -
                        points.coords.vec(simplices.vertices.vec(s, 0), i);
            }

            blift = ab[X] * ab[X] + ab[Y] * ab[Y] + ab[Z] * ab[Z];
            clift = ac[X] * ac[X] + ac[Y] * ac[Y] + ac[Z] * ac[Z];
            dlift = ad[X] * ad[X] + ad[Y] * ad[Y] + ad[Z] * ad[Z];

            simplices.subdets.store(s, 0, determinant3x3(ab[Y], ab[Z], blift,
                                                         ac[Y], ac[Z], clift,
                                                         ad[Y], ad[Z], dlift));

            simplices.subdets.store(s, 1, determinant3x3(ab[X], ab[Z], blift,
                                                         ac[X], ac[Z], clift,
                                                         ad[X], ad[Z], dlift));

            simplices.subdets.store(s, 2, determinant3x3(ab[X], ab[Y], blift,
                                                         ac[X], ac[Y], clift,
                                                         ad[X], ad[Y], dlift));

            simplices.subdets.store(s, 3, determinant3x3(ab[X], ab[Y], ab[Z],
                                                         ac[X], ac[Y], ac[Z],
                                                         ad[X], ad[Y], ad[Z]));

        }

    }

#endif

    template<class SimplexArray, class PointArray>
    static Precision
    insphere_fast(const tIndexType &s, const tIndexType &pe, const SimplexArray &simplices, const PointArray &points) {

        const tDimType D = SimplexArray::D;

        if constexpr (SimplexArray::hasSubdets) {
            Precision ae[D];
            Precision elift;

            for (uint i = 0; i < D; ++i) {
                ae[i] = points.coords(pe, i) - points.coords(simplices.vertices(s, 0), i);
            }

            elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

            return -ae[X] * simplices.subdets(s, 0) + ae[Y] * simplices.subdets(s, 1) -
                   ae[Z] * simplices.subdets(s, 2) +
                   elift * simplices.subdets(s, 3);
        } else {
            Precision ae[D], be[D], ce[D], de[D];
            Precision alift, blift, clift, dlift;
            Precision ab, bc, cd, da, ac, bd;
            Precision abc, bcd, cda, dab;

            for (uint i = 0; i < D; ++i) {
                ae[i] = points.coords(simplices.vertices(s, 0), i) - points.coords(pe, i);
                be[i] = points.coords(simplices.vertices(s, 1), i) - points.coords(pe, i);
                ce[i] = points.coords(simplices.vertices(s, 2), i) - points.coords(pe, i);
                de[i] = points.coords(simplices.vertices(s, 3), i) - points.coords(pe, i);
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

    template<class SimplexArray, class PointArray>
    static std::tuple<Precision, Precision>
    insphere(const tIndexType &s, const tIndexType &pe, const SimplexArray &simplices, const PointArray &points) {

        const tDimType D = SimplexArray::D;

        if constexpr (SimplexArray::hasSubdets) {
            Precision ae[D];
            Precision elift;
            Precision det;

            for (uint i = 0; i < D; ++i) {
                ae[i] = points.coords(pe, i) - points.coords(simplices.vertices(s, 0), i);
            }

            elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

            det = -ae[X] * simplices.subdets(s, 0) + ae[Y] * simplices.subdets(s, 1) -
                  ae[Z] * simplices.subdets(s, 2) +
                  elift * simplices.subdets(s, 3);

            return std::make_tuple(det, 0);
        } else {
            Precision ae[D], be[D], ce[D], de[D];
            Precision aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
            Precision aexcey, cexaey, bexdey, dexbey;
            Precision alift, blift, clift, dlift;
            Precision ab, bc, cd, da, ac, bd;
            Precision abc, bcd, cda, dab;
            Precision aezplus, bezplus, cezplus, dezplus;
            Precision aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
            Precision cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
            Precision aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
            Precision det;
            Precision permanent;

            for (uint i = 0; i < D; ++i) {
                ae[i] = points.coords(simplices.vertices(s, 0), i) - points.coords(pe, i);
                be[i] = points.coords(simplices.vertices(s, 1), i) - points.coords(pe, i);
                ce[i] = points.coords(simplices.vertices(s, 2), i) - points.coords(pe, i);
                de[i] = points.coords(simplices.vertices(s, 3), i) - points.coords(pe, i);
            }

            aexbey = ae[X] * be[Y];
            bexaey = be[X] * ae[Y];
            ab = aexbey - bexaey;

            bexcey = be[X] * ce[Y];
            cexbey = ce[X] * be[Y];
            bc = bexcey - cexbey;

            cexdey = ce[X] * de[Y];
            dexcey = de[X] * ce[Y];
            cd = cexdey - dexcey;

            dexaey = de[X] * ae[Y];
            aexdey = ae[X] * de[Y];
            da = dexaey - aexdey;

            aexcey = ae[X] * ce[Y];
            cexaey = ce[X] * ae[Y];
            ac = aexcey - cexaey;

            bexdey = be[X] * de[Y];
            dexbey = de[X] * be[Y];
            bd = bexdey - dexbey;

            abc = ae[Z] * bc - be[Z] * ac + ce[Z] * ab;
            bcd = be[Z] * cd - ce[Z] * bd + de[Z] * bc;
            cda = ce[Z] * da + de[Z] * ac + ae[Z] * cd;
            dab = de[Z] * ab + ae[Z] * bd + be[Z] * da;

            alift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];
            blift = be[X] * be[X] + be[Y] * be[Y] + be[Z] * be[Z];
            clift = ce[X] * ce[X] + ce[Y] * ce[Y] + ce[Z] * ce[Z];
            dlift = de[X] * de[X] + de[Y] * de[Y] + de[Z] * de[Z];

            det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

            aezplus = Absolute(ae[Z]);
            bezplus = Absolute(be[Z]);
            cezplus = Absolute(ce[Z]);
            dezplus = Absolute(de[Z]);
            aexbeyplus = Absolute(aexbey);
            bexaeyplus = Absolute(bexaey);
            bexceyplus = Absolute(bexcey);
            cexbeyplus = Absolute(cexbey);
            cexdeyplus = Absolute(cexdey);
            dexceyplus = Absolute(dexcey);
            dexaeyplus = Absolute(dexaey);
            aexdeyplus = Absolute(aexdey);
            aexceyplus = Absolute(aexcey);
            cexaeyplus = Absolute(cexaey);
            bexdeyplus = Absolute(bexdey);
            dexbeyplus = Absolute(dexbey);
            permanent = ((cexdeyplus + dexceyplus) * bezplus
                         + (dexbeyplus + bexdeyplus) * cezplus
                         + (bexceyplus + cexbeyplus) * dezplus)
                        * alift
                        + ((dexaeyplus + aexdeyplus) * cezplus
                           + (aexceyplus + cexaeyplus) * dezplus
                           + (cexdeyplus + dexceyplus) * aezplus)
                          * blift
                        + ((aexbeyplus + bexaeyplus) * dezplus
                           + (bexdeyplus + dexbeyplus) * aezplus
                           + (dexaeyplus + aexdeyplus) * bezplus)
                          * clift
                        + ((bexceyplus + cexbeyplus) * aezplus
                           + (cexaeyplus + aexceyplus) * bezplus
                           + (aexbeyplus + bexaeyplus) * cezplus)
                          * dlift;

            return std::make_tuple(det, permanent);
        }

    }

    template<class SimplexArray, class PointArray>
    static Precision
    insphere_adapt(const Precision permanent, const tIndexType &s, const tIndexType &pe, const SimplexArray &simplices,
                   const PointArray &points) {

        const tDimType D = SimplexArray::D;
        Precision a[D], b[D], c[D], d[D], e[D];

        for (uint i = 0; i < D; ++i) {
            a[i] = points.coords(simplices.vertices(s, 0), i);
            b[i] = points.coords(simplices.vertices(s, 1), i);
            c[i] = points.coords(simplices.vertices(s, 2), i);
            d[i] = points.coords(simplices.vertices(s, 3), i);
            e[i] = points.coords(pe, i);
        }

        return super::insphere_adapt(a, b, c, d, e, permanent);

    }

    template<class SimplexArray, class PointArray>
    static Precision
    pred_insphere(const tIndexType &s, const tIndexType &pe, const SimplexArray &simplices,
                  const PointArray &points) {

        const tDimType D = SimplexArray::D;
        Precision a[D], b[D], c[D], d[D], e[D];

        for (uint i = 0; i < D; ++i) {
            a[i] = points.coords(simplices.vertices(s, 0), i);
            b[i] = points.coords(simplices.vertices(s, 1), i);
            c[i] = points.coords(simplices.vertices(s, 2), i);
            d[i] = points.coords(simplices.vertices(s, 3), i);
            e[i] = points.coords(pe, i);
        }

        return super::insphere(a, b, c, d, e);

    }

    template<class SimplexArray, class PointArray>
    static Precision
    pred_orient(const tIndexType &s, const SimplexArray &simplices,
                  const PointArray &points) {

        const tDimType D = SimplexArray::D;
        Precision a[D], b[D], c[D], d[D];

        for (uint i = 0; i < D; ++i) {
            a[i] = points.coords(simplices.vertices(s, 0), i);
            b[i] = points.coords(simplices.vertices(s, 1), i);
            c[i] = points.coords(simplices.vertices(s, 2), i);
            d[i] = points.coords(simplices.vertices(s, 3), i);
        }

        return super::orient(a, b, c, d);

    }

#ifdef HAS_Vc

    template<class SimplexArray, class PointArray, typename Simplex>
    static Vc::Vector<Precision>
    insphere_fast(const Simplex &s, const Vc::Vector<tIndexType> &pe, const SimplexArray &simplices,
                  const PointArray &points) {

        const tDimType D = SimplexArray::D;

        if constexpr (SimplexArray::hasSubdets) {
            Vc::Vector<Precision> ae[D];
            Vc::Vector<Precision> elift;

            for (uint i = 0; i < D; ++i) {
                ae[i] = points.coords.vec(pe, i) - points.coords.vec(simplices.vertices.vec(s, 0), i);
            }

            elift = ae[X] * ae[X] + ae[Y] * ae[Y] + ae[Z] * ae[Z];

            return -ae[X] * simplices.subdets.vec(s, 0) + ae[Y] * simplices.subdets.vec(s, 1) -
                   ae[Z] * simplices.subdets.vec(s, 2) +
                   elift * simplices.subdets.vec(s, 3);
        } else {
            Vc::Vector<Precision> ae[D], be[D], ce[D], de[D];
            Vc::Vector<Precision> alift, blift, clift, dlift;
            Vc::Vector<Precision> ab, bc, cd, da, ac, bd;
            Vc::Vector<Precision> abc, bcd, cda, dab;

            for (uint i = 0; i < D; ++i) {
                ae[i] = points.coords.vec(simplices.vertices.vec(s, 0), i) - points.coords.vec(pe, i);
                be[i] = points.coords.vec(simplices.vertices.vec(s, 1), i) - points.coords.vec(pe, i);
                ce[i] = points.coords.vec(simplices.vertices.vec(s, 2), i) - points.coords.vec(pe, i);
                de[i] = points.coords.vec(simplices.vertices.vec(s, 3), i) - points.coords.vec(pe, i);
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
};

#endif

