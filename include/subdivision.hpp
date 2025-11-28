#pragma once

#ifndef SUBDIVISION_HPP
#define SUBDIVISION_HPP

#include "uni_b1.hpp"
#include <unistd.h>
#include <getopt.h>
#include <iomanip> 
// #include <H5Cpp.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_naive_point_location.h>

// using namespace H5;
// To-Do: Implement hdf5 output for the skyscraper invariant
using namespace graded_linalg;

namespace hnf{

struct Slope_subdivision;

using K = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = K::Point_3; // 3d point
using Point_2 = K::Point_2; // 2d point
using Segment_2 = K::Segment_2; // 2d line segment
using Polyhedron = CGAL::Polyhedron_3<K>; // 3d Polyhedron stores the lower envelope / lower convex hull
using Vector_3 = K::Vector_3; // 3d vector

// Track which segment in the arrangement came from which segment in the hull
struct Segment_w_labels {
    Segment_2 segment;
    std::optional<size_t> vertex_idx_left; 
    std::optional<size_t> vertex_idx_right;
};


struct face_data{
    int subspace_index;
    std::array<double, 3> slope_polynomial;
    std::unique_ptr<Uni_B1> submodule;
    std::unique_ptr<Uni_B1> quotient;
    
    face_data() : subspace_index(0) {}
    
    // Copy constructor
    face_data(const face_data& other) 
        : subspace_index(other.subspace_index) {
        if (other.submodule){
            submodule = std::make_unique<Uni_B1>(*other.submodule);
        } 
        if (other.quotient){
            quotient = std::make_unique<Uni_B1>(*other.quotient);
        }
    }
    
    // Move constructor
    face_data(face_data&& other) = default;
    
    // Copy assignment
    face_data& operator=(const face_data& other) {
        if (this != &other) {
            subspace_index = other.subspace_index;
            if (other.submodule){
                submodule = std::make_unique<Uni_B1>(*other.submodule);
            } 
            if (other.quotient){
                quotient = std::make_unique<Uni_B1>(*other.quotient);
            }
        }
        return *this;
    }
    
    // Move assignment
    face_data& operator=(face_data&& other) = default;
};


using Traits = CGAL::Arr_segment_traits_2<K>;


using Arrangement = CGAL::Arrangement_2<
    CGAL::Arr_segment_traits_2<K>,
    CGAL::Arr_extended_dcel<
        CGAL::Arr_segment_traits_2<K>,
        CGAL::Arr_vertex_base<typename CGAL::Arr_segment_traits_2<K>::Point_2>,
        CGAL::Arr_halfedge_base<typename CGAL::Arr_segment_traits_2<K>::X_monotone_curve_2>,
        face_data
    >
>;

struct Slope_subdivision {
    Arrangement arr;
    Slope_subdivision() = default;
    Slope_subdivision(Arrangement arrangement) : arr(std::move(arrangement)) {}

    void export_to_svg(const std::string& filename, 
        double axes_origin_x = 0.0, 
        double axes_origin_y = 0.0) const;
};

Slope_subdivision compute_slope_subdivision(Uni_B1& res, 
    const pair<r2degree>& bounds, 
    const vec<vec<SparseMatrix<int>>>& subspaces,
    const r2degree& cell_start,
    const r2degree& cell_boundary);

std::vector<Point_3> dual_points_polys(const vec<std::array<double,3>>& polynomials);


Arrangement subdivision_from_polynomials(const vec<std::array<double,3>>& polynomials);

} // namespace hnf

#endif // subdivision.hpp