#include "hnf_uni_b1.hpp"

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


void Uni_B1::compute_arrangement_quotients(vec<SparseMatrix<int>> subspaces){
    //TO-DO implement
};


void Uni_B1::compute_slope_subdivision(const pair<r2degree>& bounds, 
    const vec<vec<SparseMatrix<int>>>& subspaces,
    const r2degree& cell_start,
    const r2degree& cell_boundary){
    auto& X = this->d1;
    int k = X.get_num_rows();
    if(k == 1){
        std::cout << "Tried to compute slope subdivision for a module of dimension 1. Nothing to do." << std::endl;
        return;
    }
    vec<std::array<double, 3>> slope_polynomials;
    if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
    }
    vec<R2Mat> submodules = vec<R2Mat>();
    X.sort_compatibly();
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

    this->slope_subdiv = std::make_unique<Slope_subdivision>(
        subdivision_from_polynomials(subspaces[k-1], slope_polynomials, submodules, cell_start, cell_boundary)
    );

    this->compute_arrangement_quotients(subspaces[k-1]);
}



Uni_B1::Uni_B1(R2Mat&& d1_, bool is_minimal)
    : d1(std::move(d1_)) {
    if(is_minimal && d1.get_num_rows() == 1){
        assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
        d2 = R2Mat(d1.get_num_cols()-1, d1.get_num_cols());
        d2.data = vec< vec<int> >(d1.get_num_cols()-1);
        d2.row_degrees = d1.col_degrees;
        d2.col_degrees = vec<r2degree>(d1.get_num_cols()-1);
        r2degree last_degree = d1.col_degrees[0];
        for(int i = 1; i < d1.get_num_cols(); i++){
            r2degree join = Degree_traits<r2degree>::join(last_degree, d1.col_degrees[i]);
            d2.data[i] = {i -1, i};
            d2.col_degrees[i] = join;
            last_degree = d1.col_degrees[i];
        }
    } else {
        auto d1_copy = d1;
        d2 = d1_copy.graded_kernel();
    }
}


Uni_B1::Uni_B1(const R2Mat& d1_, bool is_minimal)
    : d1(d1_) {
    if(is_minimal && d1.get_num_rows() == 1){
        assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
        d2 = R2Mat(d1.get_num_cols()-1, d1.get_num_cols());
        d2.data = vec< vec<int> >(d1.get_num_cols()-1);
        d2.row_degrees = d1.col_degrees;
        d2.col_degrees = vec<r2degree>(d1.get_num_cols()-1);
        r2degree last_degree = d1.col_degrees[0];
        for(int i = 1; i < d1.get_num_cols(); i++){
            r2degree join = Degree_traits<r2degree>::join(last_degree, d1.col_degrees[i]);
            d2.data[i] = {i -1, i};
            d2.col_degrees[i] = join;
            last_degree = d1.col_degrees[i];
        }
    } else {
        auto d1_copy = d1;
        d2 = d1_copy.graded_kernel();
    }
}


int Uni_B1::dim_at(r2degree alpha) {
    int num_chains_0 = d1.num_rows_before(alpha);
    int num_chains_1 = d1.num_cols_before(alpha);
    int num_chains_2 = d2.num_cols_before(alpha);
    return num_chains_0 - num_chains_1 + num_chains_2;
}


double Uni_B1::area() const {
    auto [min, max] = d1.bounding_box();
    double base_area = d1.get_num_rows()*(max.first - min.first) * (max.second - min.second);
    for(const auto& degree : d1.col_degrees){
        base_area -= (max.first - degree.first) * (max.second - degree.second);
    }
    for(const auto& degree : d2.col_degrees){
        base_area += (max.first - degree.first) * (max.second - degree.second);
    }
    return base_area;
}


double Uni_B1::area(const r2degree& bound) const {
    auto [min, max] = d1.bounding_box();
    max = bound;
    double base_area = d1.get_num_rows()*(max.first - min.first) * (max.second - min.second);
    for(const auto& degree : d1.col_degrees){
        base_area -= (max.first - degree.first) * (max.second - degree.second);
    }
    for(const auto& degree : d2.col_degrees){
        base_area += (max.first - degree.first) * (max.second - degree.second);
    }
    return base_area;
}


double Uni_B1::area(const pair<r2degree>& bounds) const {
    auto [min, max] = d1.bounding_box();
    max = bounds.second;
    assert(max.first >= min.first);
    assert(max.second >= min.second);
    assert(min == d1.row_degrees[0]);
    r2degree range = bounds.second-bounds.first;
    double range_area = range.first * range.second;
    double base_area = d1.get_num_rows()*(max.first - min.first)*(max.second - min.second);
    for(const auto& degree : d1.col_degrees){
        assert(degree.first <= max.first);
        assert(degree.second <= max.second);
        assert(degree.first >= min.first);
        assert(degree.second >= min.second);
        base_area -= (max.first - degree.first) * (max.second - degree.second);
    }
    for(const auto& degree : d2.col_degrees){
        assert(degree.first <= max.first);
        assert(degree.second <= max.second);
        assert(degree.first >= min.first);
        assert(degree.second >= min.second);
        base_area += (max.first - degree.first) * (max.second - degree.second);
    }
    double normalised_area = base_area/range_area;
    return normalised_area;
}


void Uni_B1::compute_area_polynomial(const pair<r2degree>& bounds) {
    int num_rows = d1.get_num_rows();
    assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
    bool test = true;

    r2degree range = bounds.second-bounds.first;
    double range_area = range.first * range.second;

    const r2degree& gen_degree = this->d1.row_degrees[0];
    area_polynomial.fill(0.0);

    r2degree gen_vector = bounds.second - gen_degree;
    area_polynomial[0] += num_rows * gen_vector.first * gen_vector.second;
    area_polynomial[1] -= num_rows * gen_vector.second;
    area_polynomial[2] -= num_rows * gen_vector.first;

    for(const auto& degree : d1.col_degrees){
        assert( Degree_traits<r2degree>::smaller_equal(degree, bounds.second));
        assert( Degree_traits<r2degree>::smaller_equal(bounds.first, degree) );
        r2degree rel_vector = bounds.second - degree;
        area_polynomial[0] -= rel_vector.first* rel_vector.second;
    }

    for(const auto& degree :  d1.col_degrees){
        if(degree.first == gen_degree.first){
            assert(bounds.second.second - degree.second >= 0);
            area_polynomial[1] += (bounds.second.second - degree.second);
            if( area_polynomial[1] > 0){
                std::cout << "Area polynomial[1] is positive: " << area_polynomial[1] << std::endl;
                this->d1.print_graded();
                std::cout << "Bounds: " << bounds.first.first << " " << bounds.first.second << " to "
                    << bounds.second.first << " " << bounds.second.second << std::endl;
                std::cout << "Range area: " << range_area << std::endl;
                assert(false);
            }
        }
        if(degree.second == gen_degree.second){
            assert(bounds.second.first - degree.first >= 0);
            area_polynomial[2] += (bounds.second.first - degree.first);
            if( test && area_polynomial[2] > 0){
                std::cout << "Area polynomial[2] is positive: " << area_polynomial[1] << std::endl;
                this->d1.print_graded();
                std::cout << "Bounds: " << bounds.first.first << " " << bounds.first.second << " to "
                    << bounds.second.first << " " << bounds.second.second << std::endl;
                std::cout << "Range area: " << range_area << std::endl;
                assert(false);
            }
        }
    }

    for(const auto& degree : d2.col_degrees){
        assert( Degree_traits<r2degree>::smaller_equal(degree, bounds.second));
        assert( Degree_traits<r2degree>::smaller_equal(bounds.first, degree) );
        r2degree rel_vector = bounds.second - degree;
        area_polynomial[0] += rel_vector.first* rel_vector.second;
    }
    for(int i = 0; i < 3; i++){
        area_polynomial[i] /= range_area;
    }
}

double Uni_B1::evaluate_area_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    return area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + x*y;
}

double Uni_B1::evaluate_slope_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    int k = this->d1.get_num_rows();
    return k / (area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + x*y);
}


double Uni_B1::slope() const {
    double area = this->area();
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. Consider passing a bound." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}


double Uni_B1::slope(const r2degree& bound) const {
    double area = this->area(bound);
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. The bound you passed is insufficient." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}


double Uni_B1::slope(const pair<r2degree>& bounds) const {
    double area = this->area(bounds);
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. The bound you passed is insufficient." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}

} // namespace hnf
