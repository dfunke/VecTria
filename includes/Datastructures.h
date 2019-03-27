//
//  Header.h
//  delaunay
//
//  Created by Daniel Funke on 21.03.19.
//  Copyright © 2019 Daniel Funke. All rights reserved.
//

#pragma once

#include <vector>
#include <array>

using tIndexType = unsigned long;
using tDimType = unsigned short;

template<typename Precision>
using tFloatVector = std::vector<Precision>;

using tIndexVector = std::vector<tIndexType>;

template<tDimType D, typename Precision>
class PointArray {
    
public:
    std::array<tFloatVector<Precision>, D> coords;
    
public:
    
    void ensure(const tIndexType & i){
        for(tDimType d = 0; d < D+1; ++d){
            if(coords[d].size() < i){
                coords[d].resize(i);
            }
        }
    }
    
};

template <tDimType D, typename Precision>
class Point {
    
private:
    
    PointArray<D, Precision> & m_pointArray;
    tIndexType m_idx;
    
public:
    
    Point(PointArray<D, Precision> & pointArray, tIndexType & idx)
    : m_pointArray(pointArray), m_idx(idx)
    {}
    
    const Precision & operator[](const tDimType & i) const {
        return m_pointArray.coords[i][m_idx];
    }
    
    Precision & operator[](const tDimType & i) {
        return m_pointArray.coords[i][m_idx];
    }
};

template <tDimType D>
class SimplexArray {
    
public:
    std::array<tIndexVector, D+1> vertices;
    std::array<tIndexVector, D+1> neighbors;
    
public:
    
    void ensure(const tIndexType & i){
        for(tDimType d = 0; d < D+1; ++d){
            if(vertices[d].size() < i){
                vertices[d].resize(i);
            }
            
            if(neighbors[d].size() < i){
                neighbors[d].resize(i);
            }
        }
    }
};

template <tDimType D>
class Simplex {
    
private:
    
    SimplexArray<D> & m_simplexArray;
    tIndexType m_idx;
    
public:
    
    Simplex(SimplexArray<D> & simplexArray, tIndexType & idx)
    : m_simplexArray(simplexArray), m_idx(idx)
    {
        m_simplexArray.ensure(m_idx);
    }
    
    const tIndexType & vertex(const tDimType & i) const {
        return m_simplexArray.vertices[i][m_idx];
    }
    
    tIndexType & vertex(const tDimType & i) {
        return m_simplexArray.vertices[i][m_idx];
    }
    
    const tIndexType & neighbor(const tDimType & i) const {
        return m_simplexArray.neighbors[i][m_idx];
    }
    
    tIndexType & neighbor(const tDimType & i) {
        return m_simplexArray.neighbors[i][m_idx];
    }
    
};
