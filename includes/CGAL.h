#pragma once


#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/point_generators_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

//2D definitions
using Vb_2 = CGAL::Triangulation_vertex_base_with_info_2<tIndexType, K>;
using Cb_2 = CGAL::Triangulation_face_base_with_info_2<tIndexType, K>;
using Tds_2 = CGAL::Triangulation_data_structure_2<Vb_2, Cb_2>;
using DT_2 = CGAL::Delaunay_triangulation_2<K, Tds_2>;
using Point_2 = DT_2::Point;
using PointWithInfo_2 = std::pair<Point_2, tIndexType>;
using Points_2 = std::vector<PointWithInfo_2>;
using SearchTraits_2 = CGAL::Spatial_sort_traits_adapter_2<K, CGAL::First_of_pair_property_map<PointWithInfo_2>>;
//        using Fh_2 = DT_2::Face_handle;
//        using Vh_2 = DT_2::Vertex_handle;
//        using Circ_2 = CGAL::Circle_2<K>;
//        using Box_2 = CGAL::Bbox_2;
//        using Vec_2 = CGAL::Vector_2<K>;

//3D definitions
using Vb_3 = CGAL::Triangulation_vertex_base_with_info_3<tIndexType, K>;
using Cb_3 = CGAL::Triangulation_cell_base_with_info_3<tIndexType, K>;
using Tds_3 = CGAL::Triangulation_data_structure_3<Vb_3, Cb_3>;
using DT_3 = CGAL::Delaunay_triangulation_3<K, Tds_3>;
using Point_3 = DT_3::Point;
using PointWithInfo_3 = std::pair<Point_3, tIndexType>;
using Points_3 = std::vector<PointWithInfo_3>;
using SearchTraits_3 = CGAL::Spatial_sort_traits_adapter_3<K, CGAL::First_of_pair_property_map<PointWithInfo_3>>;
//        using Fh_3 = DT_3::Cell_handle;
//        using Vh_3 = DT_3::Vertex_handle;
//        using Circ_3 = CGAL::Sphere_3<K>;
//        using Box_3 = CGAL::Bbox_3;
//        using Vec_3 = CGAL::Vector_3<K>;