#pragma once

#include <random>

#include "Datastructures.h"
#include "GeometryStructures.h"
#include "CGAL.h"

#define SEED 1986

struct RandomDevice {
    static inline std::mt19937 generator = std::mt19937(SEED);
};

template<tDimType D, typename Precision>
struct Generator;

template<typename Precision>
struct Generator<2, Precision> {

    Points_2 generate(const tIndexType &n) {

        std::uniform_real_distribution<Precision> distribution;
        auto rand = std::bind(distribution, RandomDevice::generator);

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
        pa.ensure(static_cast<tIndexType>(points.size()));
        for (const auto &cp : points) {
            auto p = pa.get(cp.second);

            p[0] = cp.first.x();
            p[1] = cp.first.y();
        }
    }
    
    template<class PointArray>
    void convert(Points_2 &points, const PointArray &pa) {
        points.clear();
        points.reserve(pa.size());
        for (tIndexType i = 0; i < pa.size(); ++i) {
            auto p = pa.get(i);

            Point_2 p_cgal(p[0], p[1]);
            points.push_back(std::make_pair(p_cgal, i));
        }
    }
};

template<typename Precision>
struct Generator<3, Precision> {

    Points_3 generate(const tIndexType &n) {

        std::uniform_real_distribution<Precision> distribution;
        auto rand = std::bind(distribution, RandomDevice::generator);

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
        pa.ensure(static_cast<tIndexType>(points.size()));
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

struct UniformMover {
  
    template<class PointArray>
    void move(PointArray &points, const tVector<PointArray::D, typename PointArray::Precision> &step) {
        for (tIndexType i = 0; i < points.size(); ++i) {
            auto p = points.get(i);
            
            for(tDimType d = 0; d < PointArray::D; ++d){
                p[d] += step[d];
            }
        }
    }
    
    static std::string name() {
        return "Uniform";
    }
    
};

struct ZonedMover {
  
    template<class PointArray>
    void move(PointArray &points, const std::vector<Box<tVector<PointArray::D, typename PointArray::Precision>>> &zones,
        const std::vector<tVector<PointArray::D, typename PointArray::Precision>> &steps) {
        
            for (tIndexType i = 0; i < points.size(); ++i) {
                auto p = points.get(i);
            
                tIndexType zone = 0;
                for(; zone < zones.size(); ++zone){
                    if(zones[zone].contains(p)){
                        break;
                    }
                }
                
            for(tDimType d = 0; d < PointArray::D; ++d){
                p[d] += steps[zone][d];
            }
        }
    }
    
    static std::string name() {
        return "Zoned";
    }
    
};

struct UniformStochasticMover {
  
    template<class PointArray>
    void move(PointArray &points, const tVector<PointArray::D, typename PointArray::Precision> &step) {
        
        std::array<std::normal_distribution<typename PointArray::Precision>, PointArray::D> dists;
        for(tDimType d = 0; d < PointArray::D; ++d){
            dists[d] = std::normal_distribution<typename PointArray::Precision>(0, std::abs(step[d]));
        }
        
        for (tIndexType i = 0; i < points.size(); ++i) {
            auto p = points.get(i);
            
            for(tDimType d = 0; d < PointArray::D; ++d){
                p[d] += std::copysign(dists[d](RandomDevice::generator), step[d]);
            }
        }
    }
    
    static std::string name() {
        return "StochasticUniform";
    }
    
};

struct ZonedStochasticMover {
  
    template<class PointArray>
    void move(PointArray &points, const std::vector<Box<tVector<PointArray::D, typename PointArray::Precision>>> &zones,
        const std::vector<tVector<PointArray::D, typename PointArray::Precision>> &steps) {
        
        std::vector<std::array<std::normal_distribution<typename PointArray::Precision>, PointArray::D>> dists;
        dists.resize(zones.size());
        
        for(tIndexType i = 0; i < zones.size(); ++i){
        for(tDimType d = 0; d < PointArray::D; ++d){
            dists[i][d] = std::normal_distribution<typename PointArray::Precision>(0, std::abs(steps[i][d]));
        }
        }
        
            for (tIndexType i = 0; i < points.size(); ++i) {
                auto p = points.get(i);
            
                tIndexType zone = 0;
                for(; zone < zones.size(); ++zone){
                    if(zones[zone].contains(p)){
                        break;
                    }
                }
                
            for(tDimType d = 0; d < PointArray::D; ++d){
                p[d] += std::copysign(dists[zone][d](RandomDevice::generator), steps[zone][d]);
            }
        }
    }
    
    static std::string name() {
        return "ZonedStochastiv";
    }
    
};
