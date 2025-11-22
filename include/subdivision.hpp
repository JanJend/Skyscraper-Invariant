#pragma once

#ifndef SUBDIVISION_HPP
#define SUBDIVISION_HPP

#include "uni_b1.hpp"
#include <unistd.h>
#include <getopt.h>
#include <iomanip> 
#include <H5Cpp.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
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

    void export_to_svg(const std::string& filename, double axes_origin_x = 0.0, double axes_origin_y = 0.0) const {
        // Compute bounding box
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
        double actual_min_x = min_x;
        double actual_min_y = min_y;
        double actual_max_x = max_x;
        double actual_max_y = max_y;
        
        // Add 10% padding
        double w = max_x - min_x;
        double h = max_y - min_y;
        min_x -= w * 0.1;
        min_y -= h * 0.1;
        max_x += w * 0.1;
        max_y += h * 0.1;
        w *= 1.2;
        h *= 1.2;
        
        // Calculate aspect ratio and adjust viewBox to make more square-like
        double aspect_ratio = w / h;
        double view_w = w;
        double view_h = h;
        double view_min_x = min_x;
        double view_min_y = min_y;
        
        // Expand the narrower dimension to make aspect ratio closer to 1:1
        if (aspect_ratio > 2.0) {
            // Too wide - expand height
            double target_h = w / 2.0;
            double h_diff = target_h - h;
            view_min_y -= h_diff / 2.0;
            view_h = target_h;
        } else if (aspect_ratio < 0.5) {
            // Too tall - expand width
            double target_w = h * 2.0;
            double w_diff = target_w - w;
            view_min_x -= w_diff / 2.0;
            view_w = target_w;
        }
        
        // Fixed pixel-based sizes
        const double PIXEL_WIDTH = 800.0;
        double pixels_per_unit = PIXEL_WIDTH / view_w;
        double stroke = 2.0 / pixels_per_unit;
        double tick_size = 10.0 / pixels_per_unit;
        double font_size = 14.0 / pixels_per_unit;
        double axis_stroke = 1.5 / pixels_per_unit;

        std::ofstream svg_file(filename);
        svg_file << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='"
                << view_min_x << " " << -view_min_y - view_h << " " << view_w << " " << view_h << "'>\n";
        
        // Only apply y-flip
        svg_file << "<g transform='scale(1, -1)'>\n";

        // Draw faces with semi-transparent colors
        std::vector<std::string> colors = {"#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A",
                                            "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E2"};
        int color_idx = 0;

        for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
            if (fit->is_unbounded()) continue;
            auto ccb = fit->outer_ccb();
            auto curr = ccb;
            svg_file << "<polygon points='";
            do {
                auto pt = curr->source()->point();
                svg_file << pt.x() << "," << pt.y() << " ";
                ++curr;
            } while (curr != ccb);
            svg_file << "' fill='" << colors[color_idx % colors.size()]
                    << "' fill-opacity='0.3' stroke='none'/>\n";
            color_idx++;
        }

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
                    << "' r='" << stroke * 2 << "' fill='black'/>\n";
        }

        // Draw coordinate axes
        std::string axis_color = "#5B7C99";
        double axis_extension = w * 0.05;
        
        // X-axis
        double x_axis_start = std::max(axes_origin_x, actual_min_x);
        double x_axis_end = actual_max_x + axis_extension;
        svg_file << "<line x1='" << x_axis_start << "' y1='" << axes_origin_y
                << "' x2='" << x_axis_end << "' y2='" << axes_origin_y
                << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke << "'/>\n";

        // Y-axis
        double y_axis_start = std::max(axes_origin_y, actual_min_y);
        double y_axis_end = actual_max_y + axis_extension;
        svg_file << "<line x1='" << axes_origin_x << "' y1='" << y_axis_start
                << "' x2='" << axes_origin_x << "' y2='" << y_axis_end
                << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke << "'/>\n";

        // X-axis ticks and labels
        int num_ticks = 3;
        double x_tick_range = x_axis_end - x_axis_start;
        double x_step = x_tick_range / num_ticks;
        for (int i = 0; i <= num_ticks; ++i) {
            double x_val = x_axis_start + i * x_step;
            svg_file << "<line x1='" << x_val << "' y1='" << axes_origin_y
                    << "' x2='" << x_val << "' y2='" << (axes_origin_y + tick_size)
                    << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke * 0.7 << "'/>\n";
            
            double text_y = axes_origin_y + tick_size * 2;
            svg_file << "<text x='" << x_val << "' y='" << text_y
                    << "' font-size='" << font_size << "' font-family='Arial, sans-serif' "
                    << "text-anchor='middle' fill='" << axis_color
                    << "' transform='scale(1, -1)' transform-origin='" << x_val << " " << text_y << "'>"
                    << std::fixed << std::setprecision(3) << x_val << "</text>\n";
        }

        // Y-axis ticks and labels
        double y_tick_range = y_axis_end - y_axis_start;
        double y_step = y_tick_range / num_ticks;
        for (int i = 0; i <= num_ticks; ++i) {
            double y_val = y_axis_start + i * y_step;
            svg_file << "<line x1='" << axes_origin_x << "' y1='" << y_val
                    << "' x2='" << (axes_origin_x + tick_size) << "' y2='" << y_val
                    << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke * 0.7 << "'/>\n";
            
            double text_x = axes_origin_x + tick_size * 2;
            svg_file << "<text x='" << text_x << "' y='" << y_val
                    << "' font-size='" << font_size << "' font-family='Arial, sans-serif' "
                    << "text-anchor='start' fill='" << axis_color
                    << "' transform='scale(1, -1)' transform-origin='" << text_x << " " << y_val << "'>"
                    << std::fixed << std::setprecision(3) << y_val << "</text>\n";
        }

        svg_file << "</g>\n</svg>";
        svg_file.close();
        std::cout << "SVG exported to " << filename << std::endl;
    }
};

std::vector<Point_3> dual_points_polys(const vec<std::array<double,3>>& polynomials);


Arrangement subdivision_from_polynomials(const vec<std::array<double,3>>& polynomials);

} // namespace hnf

#endif // subdivision.hpp