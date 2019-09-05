#pragma once

#include "Datastructures.h"
#include "GeometryStructures.h"
#include "CGAL.h"

template<tDimType D>
struct Triangulator;

template<>
struct Triangulator<2> {

    auto cgal(const Points_2 &P) {


        DT_2 T(P.begin(), P.end());

        return T;
    }

    template<class SimplexArray, class DT>
    auto convert(const DT &T) {

        tIndexType simplexId = 0;
        for (auto it = T.finite_faces_begin(); it != T.finite_faces_end(); ++it) {
            it->info() = simplexId++;
        }

        SimplexArray simplices;
        simplices.ensure(simplexId - 1);

        for (auto it = T.finite_faces_begin(); it != T.finite_faces_end(); ++it) {
            tIndexType id = it->info();
            auto s = simplices.get(id);

            s.vertex(0) = it->vertex(0)->info();
            s.vertex(1) = it->vertex(1)->info();
            s.vertex(2) = it->vertex(2)->info();

            s.neighbor(0) = T.is_infinite(it->neighbor(0)) ? INF : it->neighbor(0)->info();
            s.neighbor(1) = T.is_infinite(it->neighbor(1)) ? INF : it->neighbor(1)->info();
            s.neighbor(2) = T.is_infinite(it->neighbor(2)) ? INF : it->neighbor(2)->info();

            if constexpr (SimplexArray::hasOppVertex) {
                s.oppVertex(0) = T.is_infinite(it->neighbor(0)) ? INF : it->neighbor(0)->vertex(
                        it->neighbor(0)->index(it))->info();
                s.oppVertex(1) = T.is_infinite(it->neighbor(1)) ? INF : it->neighbor(1)->vertex(
                        it->neighbor(1)->index(it))->info();
                s.oppVertex(2) = T.is_infinite(it->neighbor(2)) ? INF : it->neighbor(2)->vertex(
                        it->neighbor(2)->index(it))->info();
            }
        }

        return simplices;

    }
};

template<>
struct Triangulator<3> {

    auto cgal(const Points_3 &P) {

        DT_3 T(P.begin(), P.end());

        return T;
    }

    template<class SimplexArray, class DT>
    auto convert(const DT &T) {

        tIndexType simplexId = 0;
        for (auto it = T.finite_cells_begin(); it != T.finite_cells_end(); ++it) {
            it->info() = simplexId++;
        }

        SimplexArray simplices;
        simplices.ensure(simplexId - 1);

        for (auto it = T.finite_cells_begin(); it != T.finite_cells_end(); ++it) {
            tIndexType id = it->info();
            auto s = simplices.get(id);

            s.vertex(0) = it->vertex(0)->info();
            s.vertex(1) = it->vertex(1)->info();
            s.vertex(2) = it->vertex(2)->info();
            s.vertex(3) = it->vertex(3)->info();

            s.neighbor(0) = T.is_infinite(it->neighbor(0)) ? INF : it->neighbor(0)->info();
            s.neighbor(1) = T.is_infinite(it->neighbor(1)) ? INF : it->neighbor(1)->info();
            s.neighbor(2) = T.is_infinite(it->neighbor(2)) ? INF : it->neighbor(2)->info();
            s.neighbor(3) = T.is_infinite(it->neighbor(3)) ? INF : it->neighbor(3)->info();


            if constexpr (SimplexArray::hasOppVertex) {
                s.oppVertex(0) = T.is_infinite(it->neighbor(0)) ? INF : it->neighbor(0)->vertex(
                        it->neighbor(0)->index(it))->info();
                s.oppVertex(1) = T.is_infinite(it->neighbor(1)) ? INF : it->neighbor(1)->vertex(
                        it->neighbor(1)->index(it))->info();
                s.oppVertex(2) = T.is_infinite(it->neighbor(2)) ? INF : it->neighbor(2)->vertex(
                        it->neighbor(2)->index(it))->info();
                s.oppVertex(3) = T.is_infinite(it->neighbor(3)) ? INF : it->neighbor(3)->vertex(
                        it->neighbor(3)->index(it))->info();
            }
        }

        return simplices;

    }
};