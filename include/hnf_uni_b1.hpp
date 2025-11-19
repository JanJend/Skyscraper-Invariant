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
    std::array<double, 3> area_polynomial;
    std::unique_ptr<Slope_subdivision<index>> slope_subdiv;  

    void compute_arrangement_quotients(vec<SparseMatrix<index>> subspaces);

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
    int subspace_index;
    std::array<double, 3> slope_polynomial;
    std::unique_ptr<Uni_B1<index>> submodule;
    std::unique_ptr<Uni_B1<index>> quotient;
    
    face_data() : subspace_index(0) {}
    
    // Copy constructor
    face_data(const face_data& other) 
        : subspace_index(other.subspace_index) {
        if (other.submodule){
            submodule = std::make_unique<Uni_B1<index>>(*other.submodule);
        } 
        if (other.quotient){
            quotient = std::make_unique<Uni_B1<index>>(*other.quotient);
        }
    }
    
    // Move constructor
    face_data(face_data&& other) = default;
    
    // Copy assignment
    face_data& operator=(const face_data& other) {
        if (this != &other) {
            subspace_index = other.subspace_index;
            if (other.submodule){
                submodule = std::make_unique<Uni_B1<index>>(*other.submodule);
            } 
            if (other.quotient){
                quotient = std::make_unique<Uni_B1<index>>(*other.quotient);
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
    Slope_subdivision() = default;
    Slope_subdivision(Arrangement<index> arrangement) : arr(std::move(arrangement)) {}

    void export_to_svg(const std::string& filename) const {
        // Compute bounding box (should be unneccesary, i could just pass the box delimiters)
        double min_x = std::numeric_limits<double>::max();
        double min_y = std::numeric_limits<double>::max();
        double max_x = std::numeric_limits<double>::lowest();
        double max_y = std::numeric_limits<double>::lowest();
        
        for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
            auto pt = vit->point();
            double x = CGAL::to_double(pt.x());
            double y = CGAL::to_double(pt.y());
            min_x = std::min(min_x, x);
            min_y = std::min(min_y, y);
            max_x = std::max(max_x, x);
            max_y = std::max(max_y, y);
        }
        
        // Add 10% padding
        double w = max_x - min_x;
        double h = max_y - min_y;
        min_x -= w * 0.1; min_y -= h * 0.1;
        w *= 1.2; h *= 1.2;
        
        double scale = std::max(w, h);
        double stroke = scale * 0.002;
        double radius = scale * 0.01;
        
        

        std::ofstream svg_file(filename);
        svg_file << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='"
             << min_x << " " << min_y << " " << w << " " << h << "'>\n";
        
        // Draw edges
        for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
            auto curve = eit->curve();
            auto src = curve.source();
            auto tgt = curve.target();
            svg_file << "<line x1='" << src.x() << "' y1='" << src.y() 
                    << "' x2='" << tgt.x() << "' y2='" << tgt.y() 
                    << "' stroke='black' stroke-width='" << stroke << "'/>\n";
        }
        
        // Draw vertices
        for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
            auto pt = vit->point();
            svg_file << "<circle cx='" << pt.x() << "' cy='" << pt.y() 
                    << "' r='" << radius << "' fill='red'/>\n";
        }
        
        svg_file << "</svg>";
        svg_file.close();
        std::cout << "SVG exported to " << filename << std::endl;
    }
};

std::vector<Point_3> dual_points_polys(const vec<std::array<double,3>>& polynomials);

template<typename index>
Arrangement<index> subdivision_from_polynomials(const vec<std::array<double,3>>& polynomials);

} // namespace hnf

#endif // HNF_UNI_B1_HPP