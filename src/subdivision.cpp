#include "subdivision.hpp"

namespace hnf {

std::vector<Point_3> dual_points_polys(const vec<std::array<double,3>>& polynomials){
    vec<Point_3> dual_points; 
    for (size_t i = 0; i < polynomials.size(); i++) {
        auto& poly = polynomials[i];
        double c = poly[0];
        double a = poly[1];
        double b = poly[2];
        dual_points.emplace_back(Point_3(a, b, -c));
    }
    return dual_points;
}

struct BoundingBox {
    double x_min, x_max, y_min, y_max;
    
    bool contains(const Point_2& p) const {
        return p.x() >= x_min && p.x() <= x_max &&
               p.y() >= y_min && p.y() <= y_max;
    }
    
    K::Iso_rectangle_2 to_cgal_rect() const {
        return K::Iso_rectangle_2(Point_2(x_min, y_min), Point_2(x_max, y_max));
    }
    
    std::vector<Segment_2> boundary_segments() const {
        return {
            Segment_2(Point_2(x_min, y_min), Point_2(x_max, y_min)),
            Segment_2(Point_2(x_max, y_min), Point_2(x_max, y_max)),
            Segment_2(Point_2(x_max, y_max), Point_2(x_min, y_max)),
            Segment_2(Point_2(x_min, y_max), Point_2(x_min, y_min))
        };
    }
};

Vector_3 get_normal_of_facet(Polyhedron::Facet_const_handle fit) {
    auto he = fit->halfedge();
    Point_3 p1 = he->vertex()->point();
    Point_3 p2 = he->next()->vertex()->point();
    Point_3 p3 = he->next()->next()->vertex()->point();
    return CGAL::cross_product(p2 - p1, p3 - p1);
}


Point_2 dual_vertex_from_normal(Vector_3& normal) {
    auto a = -normal.x() / normal.z();
    auto b = -normal.y() / normal.z();
    return Point_2(a, b);
}

Point_3 dual_3_vertex_from_normal(const Vector_3& normal, const Point_3& point_on_plane) {
    auto a = -normal.x() / normal.z();
    auto b = -normal.y() / normal.z();
    auto d = normal.x() * point_on_plane.x()
           + normal.y() * point_on_plane.y()
           + normal.z() * point_on_plane.z();
    auto c = d / normal.z();
    return Point_3(a, b, -c);
}

Point_2 dual_vertex_from_facet(Polyhedron::Facet_const_handle& fit) { 
    auto normal = get_normal_of_facet(fit);
    return dual_vertex_from_normal(normal);
}

bool is_upper_normal(const Vector_3& normal) {  
    return normal.z() > 0.000001; // Maybe even increase this number
}

bool is_upper_facet(Polyhedron::Facet_const_handle fit) {  
    return is_upper_normal(get_normal_of_facet(fit)); // Maybe even increase this number
}


std::optional<Segment_2> clip_segment_to_box(const Point_2& p1, const Point_2& p2, 
                                               const BoundingBox& box) {
    auto result = CGAL::intersection(Segment_2(p1, p2), box.to_cgal_rect());
    if (result) {
        if (const Segment_2* seg = boost::get<Segment_2>(&*result)) {
            return *seg;
        }
    }
    return std::nullopt;
}


std::optional<Segment_2> clip_ray_to_box(const Point_2& source, const K::Vector_2& dir,
                                          const BoundingBox& box) {
    K::Ray_2 ray(source, dir);
    auto result = CGAL::intersection(ray, box.to_cgal_rect());
    if (result) {
        if (const Segment_2* seg = boost::get<Segment_2>(&*result)) {
            return *seg;
        }
    }
    return std::nullopt;
}


struct DualEdges {
    std::vector<std::pair<Point_2, Point_2>> finite_segments;
    std::vector<std::pair<Point_2, K::Vector_2>> rays;
};

DualEdges collect_dual_edges_and_rays(const Polyhedron& hull, const std::map<Point_3, size_t>& point_to_index, bool info = true) {
    DualEdges result;
    std::set<std::pair<Polyhedron::Vertex_const_handle, 
                       Polyhedron::Vertex_const_handle>> processed_edges;
    
    for (auto fit = hull.facets_begin(); fit != hull.facets_end(); ++fit) {
        Vector_3 normal = get_normal_of_facet(fit);
        if (!is_upper_normal(normal)) continue;
        
        Point_2 arr_vertex = dual_vertex_from_normal(normal);

        if( info ) {
            Point_3 dual_vertex_3d = dual_3_vertex_from_normal(normal, fit->halfedge()->vertex()->point());
            std::cerr << "Now processing face with dual vertex: (" << dual_vertex_3d.x() << "," << dual_vertex_3d.y() << "," << dual_vertex_3d.z() << ")" << std::endl;
        }
        auto he_circ = fit->halfedge();
        
        do {
            // Check if edge already processed
            auto v1 = he_circ->vertex();
            auto v2 = he_circ->opposite()->vertex();
            auto edge_key = std::minmax(v1, v2);
            
            if (!processed_edges.insert(edge_key).second) {
                he_circ = he_circ->next();
                continue;
            }
            
            auto opposite_facet = he_circ->opposite()->facet();
            if (opposite_facet == Polyhedron::Facet_handle()) {
                if (info) std::cerr << "Warning: Facet has no opposite facet." << std::endl;
                he_circ = he_circ->next();
                continue;
            }
            
            Vector_3 opposite_normal = get_normal_of_facet(opposite_facet);

            if( info ) {
                if(opposite_normal.z() != 0 ) {
                    Point_3 opposite_dual_vertex_3d = dual_3_vertex_from_normal(opposite_normal, he_circ->opposite()->vertex()->point());
                    std::cerr << "  Opposite face with dual vertex: (" << opposite_dual_vertex_3d.x() << "," << opposite_dual_vertex_3d.y() << "," << opposite_dual_vertex_3d.z() << ")" << std::endl;
                    } 
                else {
                    std::cout << "  Opposite face is vertical. Normal and a point on it are " <<
                                "(" << opposite_normal.x() << "," << opposite_normal.y() << "," << opposite_normal.z() << ") and (" <<
                                he_circ->opposite()->vertex()->point().x() << "," << he_circ->opposite()->vertex()->point().y() << "," <<
                                he_circ->opposite()->vertex()->point().z() << ")\n";
                }
                
            }

            if (is_upper_normal(opposite_normal)) {
                // Finite segment case
                Point_3 p1 = he_circ->vertex()->point();
                Point_3 p2 = he_circ->opposite()->vertex()->point();
                Point_2 arr_vertex_opp = dual_vertex_from_normal(opposite_normal);
                
                if (info) {
                    std::cerr << "      Hull edge: (" << p1.x() << "," << p1.y() << "," << p1.z()
                              << ") - (" << p2.x() << "," << p2.y() << "," << p2.z() << ")" << std::endl;
                    std::cerr << "      Dual segment: (" << arr_vertex.x() << "," << arr_vertex.y()
                              << ") - (" << arr_vertex_opp.x() << "," << arr_vertex_opp.y() << ")" << std::endl;
                }
                
                result.finite_segments.push_back({arr_vertex, arr_vertex_opp});
            } else {
                // Ray case
                Point_3 p1 = he_circ->opposite()->vertex()->point(); 
                Point_3 p2 = he_circ->vertex()->point();
                Point_3 p3 = he_circ->next()->vertex()->point();
                p3 = Point_3(p3.x(), p3.y(), p3.z() + 1.0);
                Vector_3 virt_opposite_normal = CGAL::cross_product(p2 - p1, p3 - p1);
                assert(virt_opposite_normal.z() > 0);
                Point_3 virt_opposite_dual_3d = dual_3_vertex_from_normal(virt_opposite_normal, p1);
                if( info ) {
                    std::cerr << "  Virtual opposite face with dual vertex: (" << virt_opposite_dual_3d.x() << "," << virt_opposite_dual_3d.y()
                                << "," << virt_opposite_dual_3d.z() << ")" << std::endl;
                }
                Point_2 virt_arr_vertex_opp = dual_vertex_from_normal(virt_opposite_normal);
                K::Vector_2 ray_direction = virt_arr_vertex_opp - arr_vertex;
                
                if (info) {
                    std::cerr << "      Hull edge (ray): (" << p1.x() << "," << p1.y() << "," << p1.z()
                              << ") - (" << p2.x() << "," << p2.y() << "," << p2.z() << ")" << std::endl;
                
                    std::cerr << "      Dual segment (ray): (" << arr_vertex.x() << "," << arr_vertex.y()
                              << ") in direction (" << ray_direction.x() << "," << ray_direction.y() << ")" << std::endl;
                
                }
                
                result.rays.push_back({arr_vertex, ray_direction});
            }
            
            he_circ = he_circ->next();
        } while (he_circ != fit->halfedge());
    }
    
    return result;
}

Segment_2 ensure_lexicographic_order(const Segment_2& s) {
    if (CGAL::compare_xy(s.source(), s.target()) == CGAL::LARGER) {
        return Segment_2(s.target(), s.source());
    }
    return s;
}

std::vector<Segment_2> clip_and_order_segments(const DualEdges& edges, 
                                                 const BoundingBox& box,
                                                 bool info = true) {
    std::vector<Segment_2> clipped_segments;
    
    // Clip finite segments
    for (const auto& [p1, p2] : edges.finite_segments) {
        if (auto seg = clip_segment_to_box(p1, p2, box)) {
            auto s = *seg;
            auto original = s;
            s = ensure_lexicographic_order(s);
            
            if (info && original != s) {
                std::cout << "Segment was " << original.source() << " to " << original.target() << "\n";
                std::cout << "Reversing direction.\n";
            }
            
            if (info) {
                std::cerr << "  Adding clipped segment: (" << s.source().x() << "," << s.source().y()
                          << ") → (" << s.target().x() << "," << s.target().y() << ")" << std::endl;
            }
            
            clipped_segments.push_back(s);
        }
    }
    
    // Clip rays
    for (const auto& [origin, direction] : edges.rays) {
        if (auto seg = clip_ray_to_box(origin, direction, box)) {
            auto s = *seg;
            auto original = s;
            s = ensure_lexicographic_order(s);
            
            if (info && original != s) {
                std::cout << "Ray segment was " << original.source() << " to " << original.target() << "\n";
                std::cout << "Reversing direction.\n";
            }
            
            if (info) {
                std::cerr << "  Adding ray (?) after clipping : (" << s.source().x() << "," << s.source().y()
                          << ") → (" << s.target().x() << "," << s.target().y() << ")" << std::endl;
            }
            
            clipped_segments.push_back(s);
        }
    }

     // Add bounding box boundary segments
    auto boundary_segs = box.boundary_segments();
    for (auto& seg : boundary_segs) {
        seg = ensure_lexicographic_order(seg);
        if (info) {
            std::cerr << "  Adding bounding box segment: (" << seg.source().x() << "," << seg.source().y()
                      << ") → (" << seg.target().x() << "," << seg.target().y() << ")" << std::endl;
        }
        clipped_segments.push_back(seg);
    }
    
    return clipped_segments;
}

std::vector<Segment_2> extract_segments_from_hull(const Polyhedron& hull,
                                                    const BoundingBox& box,
                                                    const std::map<Point_3, size_t>& point_to_index,
                                                    bool info = true) {
    auto edges = collect_dual_edges_and_rays(hull, point_to_index, info);
    return clip_and_order_segments(edges, box, info);
}


void assign_face_data(Arrangement& arr, const Polyhedron& hull,
                        const vec<SparseMatrix<int>>& subspaces,
                     const std::vector<std::array<double,3>>& polynomials,
                     const vec<R2Mat>& submodules,
                     const std::map<Point_3, size_t>& point_to_index) {
    CGAL::Arr_naive_point_location<Arrangement> pl(arr);
    
    //To-Do
}


Arrangement subdivision_from_polynomials( const vec<SparseMatrix<int>>& subspaces,
                                                const vec<std::array<double,3>>& polynomials,
                                                 const vec<R2Mat>& submodules,
                                                 const r2degree& cell_start, 
                                                 const r2degree& cell_end, const bool info = true) {
    BoundingBox box{cell_start.first, cell_end.first, 
                           cell_start.second, cell_end.second};
    // Build convex hull in dual space
    std::vector<Point_3> dual_points = dual_points_polys(polynomials);
    std::map<Point_3, size_t> point_to_index;
    for (size_t i = 0; i < dual_points.size(); ++i) {
        point_to_index[dual_points[i]] = i;
    }
    
    Polyhedron hull;
    CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), hull);
    
    // Extract and clip segments from hull
    std::vector<Segment_2> clipped_segments = extract_segments_from_hull(hull, box, point_to_index, info);

    if(info){
        for (const auto& seg : clipped_segments) {
            std::cout << "Segment: (" << seg.source() << ") -> (" << seg.target() << ")\n";
        }
    }
    // Build arrangement with bounding box
    Arrangement arr;
    auto bbox_segments = box.boundary_segments();
    if(info){
        for (const auto& seg : bbox_segments) {
            std::cout << "Segment: (" << seg.source() << ") -> (" << seg.target() << ")\n";
        }
    }

    CGAL::insert(arr, bbox_segments.begin(), bbox_segments.end());
    CGAL::insert(arr, clipped_segments.begin(), clipped_segments.end());
    
    // Assign face data
    assign_face_data(arr, hull, subspaces, polynomials, submodules, point_to_index);
    
    return arr;
}


void compute_arrangement_quotients(vec<SparseMatrix<int>> subspaces){
    //TO-DO implement
};


Slope_subdivision compute_slope_subdivision(Uni_B1& res, 
    const pair<r2degree>& bounds, 
    const vec<vec<SparseMatrix<int>>>& subspaces,
    const r2degree& cell_start,
    const r2degree& cell_boundary){
    auto& X = res.d1;
    int k = X.get_num_rows();
    if(k == 1){
        std::cout << "Tried to compute slope subdivision for a module of dimension 1. Nothing to do." << std::endl;
        return Slope_subdivision(Arrangement());
    }
    vec<std::array<double, 3>> slope_polynomials;
    if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
    }
    vec<R2Mat> submodules = vec<R2Mat>();
    for(size_t i = 1; i < subspaces[k-1].size(); i++){
        //skip i = 0, because it is the empty space
        auto ungraded_subspace = subspaces[k-1][i];
        int num_gens = ungraded_subspace.get_num_cols();
        R2Mat subspace = R2Mat(ungraded_subspace);
        subspace.row_degrees = X.row_degrees;
        subspace.col_degrees = vec<r2degree>(num_gens, X.row_degrees[0]);
            assert(subspace.get_num_rows() == X.get_num_rows());
            assert(subspace.get_num_cols() == num_gens);
        submodules.emplace_back(X.submodule_generated_by(subspace));
        X.column_reduction_graded(); //full minimisation should not be necessary
        Uni_B1 res(submodules.back());
        res.compute_area_polynomial(bounds);  // Compute first
        slope_polynomials.emplace_back(res.area_polynomial);
        for(auto& coeff : slope_polynomials.back()){
            coeff /= static_cast<double>(num_gens);
        }
    }

    Arrangement arr = subdivision_from_polynomials(subspaces[k-1], slope_polynomials, submodules, cell_start, cell_boundary);
    Slope_subdivision result(arr);
    compute_arrangement_quotients(subspaces[k-1]);

    return result;
}

void Slope_subdivision::export_to_svg(const std::string& filename,
    double axes_origin_x,
    double axes_origin_y) const {
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
    
    // Calculate aspect ratio and determine scaling factors
    double aspect_ratio = w / h;
    double scale_x = 1.0;
    double scale_y = 1.0;
    
    if (aspect_ratio > 2.0) {
        // Too wide - expand y
        scale_y = aspect_ratio / 2.0;
    } else if (aspect_ratio < 0.5) {
        // Too tall - expand x
        scale_x = aspect_ratio * 2.0;
    }
    
    // Apply scaling to dimensions and bounds
    double scaled_w = w * scale_x;
    double scaled_h = h * scale_y;
    double scaled_min_x = min_x * scale_x;
    double scaled_min_y = min_y * scale_y;
    double scaled_max_x = max_x * scale_x;
    double scaled_max_y = max_y * scale_y;
    double scaled_actual_min_x = actual_min_x * scale_x;
    double scaled_actual_max_x = actual_max_x * scale_x;
    double scaled_actual_min_y = actual_min_y * scale_y;
    double scaled_actual_max_y = actual_max_y * scale_y;
    double scaled_axes_origin_x = axes_origin_x * scale_x;
    double scaled_axes_origin_y = axes_origin_y * scale_y;
    
    // Fixed pixel-based sizes
    const double PIXEL_WIDTH = 800.0;
    double pixels_per_unit = PIXEL_WIDTH / scaled_w;
    double stroke = 2.0 / pixels_per_unit;
    double tick_size = 10.0 / pixels_per_unit;
    double font_size = 14.0 / pixels_per_unit;
    double axis_stroke = 1.5 / pixels_per_unit;

    std::ofstream svg_file(filename);
    svg_file << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='"
            << scaled_min_x << " " << -scaled_min_y - scaled_h << " " << scaled_w << " " << scaled_h << "'>\n";
    // Only apply y-flip
    svg_file << "<g transform='scale(1, -1)'>\n";

    // Draw faces with semi-transparent colors
    std::vector<std::string> colors = {
        "#FF6B6B", "#4ECDC4", "#45B7D1", "#FFA07A", "#98D8C8", "#F7DC6F", "#BB8FCE", "#85C1E2",
        "#F8B195", "#F67280", "#C06C84", "#6C5B7B", "#355C7D", "#99B898", "#FECEA8", "#FF847C",
        "#E84A5F", "#2A363B", "#A8E6CF", "#DCEDC1", "#FFD3B6", "#FFAAA5", "#AA96DA", "#FCBAD3",
        "#FFFFD2", "#A8DADC", "#457B9D", "#1D3557", "#E63946", "#F1FAEE", "#A8E6CE", "#FFE156"
    };
    int color_idx = 0;

    for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
        if (fit->is_unbounded()) continue;
        auto ccb = fit->outer_ccb();
        auto curr = ccb;
        svg_file << "<polygon points='";
        do {
            auto pt = curr->source()->point();
            double x = CGAL::to_double(pt.x()) * scale_x;
            double y = CGAL::to_double(pt.y()) * scale_y;
            svg_file << x << "," << y << " ";
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
        double x1 = CGAL::to_double(src.x()) * scale_x;
        double y1 = CGAL::to_double(src.y()) * scale_y;
        double x2 = CGAL::to_double(tgt.x()) * scale_x;
        double y2 = CGAL::to_double(tgt.y()) * scale_y;
        svg_file << "<line x1='" << x1 << "' y1='" << y1
                << "' x2='" << x2 << "' y2='" << y2
                << "' stroke='black' stroke-width='" << stroke << "'/>\n";
    }
    
    // Draw vertices
    for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
        auto pt = vit->point();
        double x = CGAL::to_double(pt.x()) * scale_x;
        double y = CGAL::to_double(pt.y()) * scale_y;
        svg_file << "<circle cx='" << x << "' cy='" << y
                << "' r='" << stroke * 2 << "' fill='black'/>\n";
    }

    // Draw coordinate axes
    std::string axis_color = "#5B7C99";
    double axis_extension = scaled_w * 0.05;
    
    // X-axis
    double x_axis_start = std::max(scaled_axes_origin_x, scaled_actual_min_x);
    double x_axis_end = scaled_actual_max_x + axis_extension;
    svg_file << "<line x1='" << x_axis_start << "' y1='" << scaled_axes_origin_y
            << "' x2='" << x_axis_end << "' y2='" << scaled_axes_origin_y
            << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke << "'/>\n";

    // Y-axis
    double y_axis_start = std::max(scaled_axes_origin_y, scaled_actual_min_y);
    double y_axis_end = scaled_actual_max_y + axis_extension;
    svg_file << "<line x1='" << scaled_axes_origin_x << "' y1='" << y_axis_start
            << "' x2='" << scaled_axes_origin_x << "' y2='" << y_axis_end
            << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke << "'/>\n";

    // X-axis ticks and labels (labels show original unscaled coordinates)
    int num_ticks = 3;
    double x_tick_range = x_axis_end - x_axis_start;
    double x_step = x_tick_range / num_ticks;
    double unscaled_x_start = x_axis_start / scale_x;
    double unscaled_x_step = x_step / scale_x;
    for (int i = 0; i <= num_ticks; ++i) {
        double x_val_scaled = x_axis_start + i * x_step;
        double x_val_unscaled = unscaled_x_start + i * unscaled_x_step;
        svg_file << "<line x1='" << x_val_scaled << "' y1='" << (scaled_axes_origin_y - tick_size)
                << "' x2='" << x_val_scaled << "' y2='" << (scaled_axes_origin_y + tick_size)
                << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke << "'/>\n";
        double text_y = scaled_axes_origin_y - tick_size * 2.5;  // Changed: minus instead of plus
        svg_file << "<text x='" << x_val_scaled << "' y='" << text_y
                << "' font-size='" << font_size << "' font-family='Arial, sans-serif' "
                << "text-anchor='middle' fill='" << axis_color
                << "' transform='scale(1, -1)' transform-origin='" << x_val_scaled << " " << text_y << "'>"
                << std::fixed << std::setprecision(2) << x_val_unscaled << "</text>\n";
    }

    // Y-axis ticks and labels (labels show original unscaled coordinates)
    double y_tick_range = y_axis_end - y_axis_start;
    double y_step = y_tick_range / num_ticks;
    double unscaled_y_start = y_axis_start / scale_y;
    double unscaled_y_step = y_step / scale_y;
    for (int i = 0; i <= num_ticks; ++i) {
        double y_val_scaled = y_axis_start + i * y_step;
        double y_val_unscaled = unscaled_y_start + i * unscaled_y_step;
        svg_file << "<line x1='" << (scaled_axes_origin_x - tick_size) << "' y1='" << y_val_scaled
                << "' x2='" << (scaled_axes_origin_x + tick_size) << "' y2='" << y_val_scaled
                << "' stroke='" << axis_color << "' stroke-width='" << axis_stroke << "'/>\n";
        double text_x = scaled_axes_origin_x - tick_size * 2.5;  // Changed: minus instead of plus
        svg_file << "<text x='" << text_x << "' y='" << y_val_scaled
                << "' font-size='" << font_size << "' font-family='Arial, sans-serif' "
                << "text-anchor='end' fill='" << axis_color  // Changed: 'end' instead of 'start'
                << "' transform='scale(1, -1)' transform-origin='" << text_x << " " << y_val_scaled << "'>"
                << std::fixed << std::setprecision(2) << y_val_unscaled << "</text>\n";
    }

    svg_file << "</g>\n</svg>";
    svg_file.close();
    std::cout << "SVG exported to " << filename << std::endl;
}

} // namespace hnf