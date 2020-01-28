#pragma once

#include <random>

#include "Datastructures.h"
#include "GeometryStructures.h"
#include "CGAL.h"

#define SEED 1986

template<tDimType D, typename Precision>
struct Generator;

template<typename Precision>
struct Generator<2, Precision> {

    Points_2 generate(const tIndexType &n) {

        std::mt19937 generator(SEED);
        std::uniform_real_distribution<Precision> distribution;
        auto rand = std::bind(distribution, generator);

        Points_2 points;
        points.reserve(n);

        for (tIndexType i = 0; i < n; ++i) {
            Point_2 p_cgal(rand(), rand());
            points.push_back(std::make_pair(p_cgal, i));
        }

        SearchTraits_2 traits;
        CGAL::hilbert_sort(points.begin(), points.end(), traits);

        return points;
    }

    template<class PointArray>
    void convert(PointArray &pa, const Points_2 &points) {
        pa.ensure(points.size());
        for (const auto &cp : points) {
            auto p = pa.get(cp.second);

            p[0] = cp.first.x();
            p[1] = cp.first.y();
        }
    }
};

template<typename Precision>
struct Generator<3, Precision> {

    Points_3 generate(const tIndexType &n) {

        std::mt19937 generator(SEED);
        std::uniform_real_distribution<Precision> distribution;
        auto rand = std::bind(distribution, generator);

        Points_3 points;
        points.reserve(n);

        for (tIndexType i = 0; i < n; ++i) {
            Point_3 p_cgal(rand(), rand(), rand());
            points.push_back(std::make_pair(p_cgal, i));
        }

        SearchTraits_3 traits;
        CGAL::hilbert_sort(points.begin(), points.end(), traits);

        return points;
    }

    template<class PointArray>
    void convert(PointArray &pa, const Points_3 &points) {
        pa.ensure(points.size());
        for (const auto &cp : points) {
            auto p = pa.get(cp.second);

            p[0] = cp.first.x();
            p[1] = cp.first.y();
            p[2] = cp.first.z();
        }
    }

    template<class PointArray>
    void convert(Points_3 &points, const PointArray &pa) {
        points.clear();
        points.reserve(pa.size());
        for (tIndexType i = 0; i < pa.size(); ++i) {
            auto p = pa.get(i);

            Point_3 p_cgal(p[0], p[1], p[2]);
            points.push_back(std::make_pair(p_cgal, i));
        }
    }
};