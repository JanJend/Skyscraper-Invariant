#pragma once

#ifndef HNF_UNI_B1_HPP
#define HNF_UNI_B1_HPP

#include "grlina/graded_linalg.hpp"
#include <unistd.h>
#include <getopt.h>
#include <H5Cpp.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_naive_point_location.h>

using namespace H5;
// To-Do: Implement hdf5 output for the skyscraper invariant
using namespace graded_linalg;

namespace hnf{

using R2Mat = R2GradedSparseMatrix<int>;

// Forward declarations
template<typename index>
struct Uni_B1;

template<typename index>
struct Slope_subdivision;


template<typename index>
struct Uni_B1{
    R2GradedSparseMatrix<index> d1;
    R2GradedSparseMatrix<index> d2;
    std::array<double, 4> area_polynomial;
    std::unique_ptr<Slope_subdivision<index>> slope_subdiv;  

    Uni_B1() = default;
    Uni_B1(R2GradedSparseMatrix<index>&& d1_, bool is_minimal = false);
    Uni_B1(const R2GradedSparseMatrix<index>& d1_, bool is_minimal = false);

    // Copy constructor
    Uni_B1(const Uni_B1& other) 
        : d1(other.d1), d2(other.d2), area_polynomial(other.area_polynomial) {
        if (other.slope_subdiv) {
            slope_subdiv = std::make_unique<Slope_subdivision<index>>(*other.slope_subdiv);
        }
    }
    
    // Move constructor
    Uni_B1(Uni_B1&& other) = default;
    
    // Copy assignment
    Uni_B1& operator=(const Uni_B1& other) {
        if (this != &other) {
            d1 = other.d1;
            d2 = other.d2;
            area_polynomial = other.area_polynomial;
            if (other.slope_subdiv) {
                slope_subdiv = std::make_unique<Slope_subdivision<index>>(*other.slope_subdiv);
            } else {
                slope_subdiv.reset();
            }
        }
        return *this;
    }
    
    // Move assignment
    Uni_B1& operator=(Uni_B1&& other) = default;

    index dim_at(r2degree alpha) ;
    double area() const;
    double area(const r2degree& bound) const;
    double area(const pair<r2degree>& bounds) const;
    void compute_area_polynomial(const pair<r2degree>& bounds);
    double evaluate_area_polynomial(r2degree d);
    double evaluate_slope_polynomial(r2degree d);
    double slope() const;
    double slope(const r2degree& bound) const;
    double slope(const pair<r2degree>& bounds) const;

    void compute_slope_subdivision(const pair<r2degree>& bounds, const vec<vec<SparseMatrix<int>>>& subspaces, const r2degree& cell_boundary);
};


using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point_3 = K::Point_3; // 3d point
using Point_2 = K::Point_2; // 2d point
using Segment_2 = K::Segment_2; // 2d line segment
using Polyhedron = CGAL::Polyhedron_3<K>; // 3d Polyhedron stores the lower envelope / lower convex hull


template<typename index>
struct face_data{
    index subspace_index;
    std::unique_ptr<Uni_B1<index>> summand;
    
    face_data() : subspace_index(0) {}
    
    // Copy constructor
    face_data(const face_data& other) 
        : subspace_index(other.subspace_index) {
        if (other.summand) {
            summand = std::make_unique<Uni_B1<index>>(*other.summand);
        }
    }
    
    // Move constructor
    face_data(face_data&& other) = default;
    
    // Copy assignment
    face_data& operator=(const face_data& other) {
        if (this != &other) {
            subspace_index = other.subspace_index;
            if (other.summand) {
                summand = std::make_unique<Uni_B1<index>>(*other.summand);
            } else {
                summand.reset();
            }
        }
        return *this;
    }
    
    // Move assignment
    face_data& operator=(face_data&& other) = default;
};

template<typename index>
using Traits = CGAL::Arr_segment_traits_2<K>;


template<typename index>
using Arrangement = CGAL::Arrangement_2<
    CGAL::Arr_segment_traits_2<K>,
    CGAL::Arr_extended_dcel<
        CGAL::Arr_segment_traits_2<K>,
        CGAL::Arr_vertex_base<typename CGAL::Arr_segment_traits_2<K>::Point_2>,
        CGAL::Arr_halfedge_base<typename CGAL::Arr_segment_traits_2<K>::X_monotone_curve_2>,
        face_data<index>
    >
>;

template<typename index>
struct Slope_subdivision {
    Arrangement<index> arr;
};

template<typename index>
Arrangement<index> subdivision_from_polynomials(const vec<std::array<double,4>>& polynomials);

} // namespace hnf

#endif // HNF_UNI_B1_HPP