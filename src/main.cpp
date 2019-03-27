//
// Created by funke on 11/23/16.
//

#include <random>
#include <functional>

#include "Datastructures.h"
#include "Predicates.h"

#define SEED 1986

template<tDimType D, typename Precision, class PointArray>
void generatePoints(PointArray &pa, const tIndexType &n) {

    std::mt19937 generator(SEED);
    std::uniform_real_distribution<Precision> distribution;
    auto rand = std::bind(distribution, generator);

    pa.ensure(n);
    for (tIndexType i = 0; i < n; ++i) {
        auto p = pa.get(i);

        for (tDimType d = 0; d < D; ++d) {
            p[d] = rand();
        }
    }
}

#define D 3
#define Precision double
#define N 1e4

int main() {

    PointAoA<D, Precision> points_aoa;
    generatePoints<D, Precision>(points_aoa, N);

    PointPA<D, Precision> points_pa;
    generatePoints<D, Precision>(points_pa, N);

    SimplexAoA<D> simplices_aoa;
    SimplexPA<D> simplices_pa;

    __attribute__((unused)) auto s_aoa = simplices_aoa.get(0);
    __attribute__((unused)) auto s_pa = simplices_pa.get(0);

    return 0;

}
