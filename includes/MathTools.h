#pragma once

template<typename Precision>
Precision determinant3x3(const Precision a00, const Precision a01, const Precision a02,
                         const Precision a10, const Precision a11, const Precision a12,
                         const Precision a20, const Precision a21, const Precision a22) {
    return (a00 * (a11 * a22 - a12 * a21) -
            a01 * (a10 * a22 - a12 * a20) +
            a02 * (a10 * a21 - a11 * a20));
}