#include "hnf_uni_b1.hpp"

namespace hnf {

using R2Mat = R2GradedSparseMatrix<int>;

struct r3_point_w_number{
    std::array<double, 3> point;
    int number;
};

std::vector<Point_3> dual_points_polys(const vec<std::array<double,4>>& polynomials){
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

template<typename index>
Arrangement<index> subdivision_from_polynomials(const vec<std::array<double,4>>& polynomials,
    const r2degree& cell_start, const r2degree& cell_end) {
    
    double x_min = cell_start.first, x_max = cell_end.first;
    double y_min = cell_start.second, y_max = cell_end.second;

    // Helper to check if point is in box
    auto in_box = [&](const Point_2& p) {
        return p.x() >= x_min && p.x() <= x_max && 
            p.y() >= y_min && p.y() <= y_max;
    };
    // Helper to clip segment to box
    auto clip_segment = [&](const Point_2& p1, const Point_2& p2) 
        -> std::optional<Segment_2> {
        // Cohen-Sutherland or simpler: use CGAL's intersection
        K::Iso_rectangle_2 box(Point_2(x_min, y_min), Point_2(x_max, y_max));
        auto result = CGAL::intersection(Segment_2(p1, p2), box);
        
        if (result) {
            if (const Segment_2* seg = boost::get<Segment_2>(&*result)) {
                return *seg;
            } else if (const Point_2* pt = boost::get<Point_2>(&*result)) {
                // Degenerate case: segment touches box at single point
                return std::nullopt;
            }
        }
        return std::nullopt;
    };

    // Helper to clip ray to box
    auto clip_ray = [&](const Point_2& source, const K::Vector_2& dir) 
        -> std::optional<Segment_2> {
        // Find intersection of ray with box
        K::Ray_2 ray(source, dir);
        K::Iso_rectangle_2 box(Point_2(x_min, y_min), Point_2(x_max, y_max));
        auto result = CGAL::intersection(ray, box);
        
        if (result) {
            if (const Segment_2* seg = boost::get<Segment_2>(&*result)) {
                return *seg;
            } else if (const Point_2* pt = boost::get<Point_2>(&*result)) {
                // Ray starts at boundary
                return std::nullopt;
            }
        }
        return std::nullopt;
    };


    std::vector<Point_3> dual_points = dual_points_polys(polynomials);
  
    std::map<Point_3, size_t> point_to_index;
    for (size_t i = 0; i < dual_points.size(); ++i) {
        point_to_index[dual_points[i]] = i;
    }
    
    Polyhedron hull;
    CGAL::convex_hull_3(dual_points.begin(), dual_points.end(), hull);

    std::vector<Segment_2> clipped_segments;
    
    for (auto fit = hull.facets_begin(); fit != hull.facets_end(); ++fit) {
        auto he = fit->halfedge();
        Point_3 p1 = he->vertex()->point();
        Point_3 p2 = he->next()->vertex()->point();
        Point_3 p3 = he->next()->next()->vertex()->point();
        
        auto normal = CGAL::cross_product(p2 - p1, p3 - p1);
        
        if (normal.z() < 0) {
            double a = -normal.x() / normal.z();
            double b = -normal.y() / normal.z();
            double c = p1.z() - a * p1.x() - b * p1.y();
            
            Point_2 arr_vertex(a, b);
            
            // Only process if vertex is in or near the box

                
            auto he_circ = fit->halfedge();
            do {
                auto opposite_facet = he_circ->opposite()->facet();
                
                if (opposite_facet != Polyhedron::Facet_handle()) {
                    auto he_opp = opposite_facet->halfedge();
                    Point_3 op1 = he_opp->vertex()->point();
                    Point_3 op2 = he_opp->next()->vertex()->point();
                    Point_3 op3 = he_opp->next()->next()->vertex()->point();
                    auto normal_opp = CGAL::cross_product(op2 - op1, op3 - op1);
                    
                    if (normal_opp.z() < 0) {
                        double a_opp = -normal_opp.x() / normal_opp.z();
                        double b_opp = -normal_opp.y() / normal_opp.z();
                        Point_2 arr_vertex_opp(a_opp, b_opp);
                        
                        // Clip segment to box
                        if (auto seg = clip_segment(arr_vertex, arr_vertex_opp)) {
                            clipped_segments.push_back(*seg);
                        }
                    } else {
                        // Ray case
                        Point_3 edge_start = he_circ->vertex()->point();
                        Point_3 edge_end = he_circ->next()->vertex()->point();
                        K::Vector_2 edge_dir(edge_end.x() - edge_start.x(), 
                                            edge_end.y() - edge_start.y());
                        K::Vector_2 ray_dir(-edge_dir.y(), edge_dir.x());
                        
                        if (normal_opp.z() > 0) {
                            ray_dir = -ray_dir;
                        }
                        
                        // Clip ray to box
                        if (auto seg = clip_ray(arr_vertex, ray_dir)) {
                            clipped_segments.push_back(*seg);
                        }
                    }
                } else {
                    std::cerr << "Warning: Facet has no opposite facet." << std::endl;
                }
                he_circ = he_circ->next();
            } while (he_circ != fit->halfedge());
            
        }
    }
    
    
    // Build arrangement with bounding box
    Arrangement<index> arr;
    std::vector<Segment_2> bbox_segments = {
        Segment_2(Point_2(x_min, y_min), Point_2(x_max, y_min)),
        Segment_2(Point_2(x_max, y_min), Point_2(x_max, y_max)),
        Segment_2(Point_2(x_max, y_max), Point_2(x_min, y_max)),
        Segment_2(Point_2(x_min, y_max), Point_2(x_min, y_min))
    };
    CGAL::insert(arr, bbox_segments.begin(), bbox_segments.end());
    CGAL::insert(arr, clipped_segments.begin(), clipped_segments.end());
    
    // Assign face data: each face corresponds to a vertex on the lower hull
    for (auto vit = hull.vertices_begin(); vit != hull.vertices_end(); ++vit) {
        Point_3 hull_vertex = vit->point();
        
        // Check if this vertex is on the lower hull (has at least one incident lower facet)
        bool is_lower = false;
        auto he_circ = vit->halfedge();
        do {
            if (he_circ->facet() != Polyhedron::Facet_handle()) {
                auto he = he_circ->facet()->halfedge();
                Point_3 p1 = he->vertex()->point();
                Point_3 p2 = he->next()->vertex()->point();
                Point_3 p3 = he->next()->next()->vertex()->point();
                auto normal = CGAL::cross_product(p2 - p1, p3 - p1);
                if (normal.z() < 0) {
                    is_lower = true;
                    break;
                }
            }
            he_circ = he_circ->next()->opposite();
        } while (he_circ != vit->halfedge());
        
        if (is_lower) {
            // Find the face in the arrangement containing point (hull_vertex.x(), hull_vertex.y())
            Point_2 query(hull_vertex.x(), hull_vertex.y());
            // Use a point location strategy
            CGAL::Arr_naive_point_location<Arrangement<index>> pl(arr);
            auto obj = pl.locate(query);
            
            if (auto* f = boost::get<Arrangement<index>::Face_handle>(&obj)) {
                // Get original index
                size_t orig_idx = point_to_index[hull_vertex];
                
                // Set face data
                (*f)->data().subspace_index = orig_idx;
                (*f)->data().slope_polynomial = polynomials[orig_idx];
                // Set other face_data fields as needed
            }
        }
    }
    
    return arr;
}

template<typename index>
void Uni_B1<index>::compute_slope_subdivision(const pair<r2degree>& bounds, 
    const vec<vec<SparseMatrix<int>>>& subspaces,
    const r2degree& cell_boundary){
    auto& X = this->d1;
    int k = X.get_num_rows();
    if(k == 1){
        std::cout << "Tried to compute slope subdivision for a module of dimension 1. Nothing to do." << std::endl;
        return;
    }
    vec<std::array<double, 4>> slope_polynomials;
    if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
    }

    for(auto ungraded_subspace : subspaces[k-1]){
        int num_gens = ungraded_subspace.get_num_cols();
        R2Mat subspace = R2Mat(ungraded_subspace);
        subspace.row_degrees = X.row_degrees;
        subspace.col_degrees = vec<r2degree>(num_gens, X.row_degrees[0]);
            assert(subspace.get_num_rows() == X.get_num_rows());
            assert(subspace.get_num_cols() == num_gens);
        R2Mat submodule_pres = X.submodule_generated_by(subspace);
        int local_dim = submodule_pres.get_num_rows();
        Uni_B1<int> res(submodule_pres);
        slope_polynomials.emplace_back( res.area_polynomial(bounds) );
        for(auto& coeff : slope_polynomials.back()){
            coeff /= static_cast<double>(local_dim);
        }
    }

    this->slope_subdivision = subdivision_from_polynomials<index>(slope_polynomials, X.row_degrees[0], X.row_degrees[0], cell_boundary);
}


template<typename index>
Uni_B1<index>::Uni_B1(R2GradedSparseMatrix<index>&& d1_, bool is_minimal)
    : d1(std::move(d1_)) {
    if(is_minimal && d1.get_num_rows() == 1){
        assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
        d2 = R2GradedSparseMatrix<index>(d1.get_num_cols()-1, d1.get_num_cols());
        d2.data = vec< vec<index> >(d1.get_num_cols()-1);
        d2.row_degrees = d1.col_degrees;
        d2.col_degrees = vec<r2degree>(d1.get_num_cols()-1);
        r2degree last_degree = d1.col_degrees[0];
        for(index i = 1; i < d1.get_num_cols(); i++){
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

template<typename index>
Uni_B1<index>::Uni_B1(const R2GradedSparseMatrix<index>& d1_, bool is_minimal)
    : d1(d1_) {
    if(is_minimal && d1.get_num_rows() == 1){
        assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
        d2 = R2GradedSparseMatrix<index>(d1.get_num_cols()-1, d1.get_num_cols());
        d2.data = vec< vec<index> >(d1.get_num_cols()-1);
        d2.row_degrees = d1.col_degrees;
        d2.col_degrees = vec<r2degree>(d1.get_num_cols()-1);
        r2degree last_degree = d1.col_degrees[0];
        for(index i = 1; i < d1.get_num_cols(); i++){
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

template<typename index>
index Uni_B1<index>::dim_at(r2degree alpha) {
    index num_chains_0 = d1.num_rows_before(alpha);
    index num_chains_1 = d1.num_cols_before(alpha);
    index num_chains_2 = d2.num_cols_before(alpha);
    return num_chains_0 - num_chains_1 + num_chains_2;
}

template<typename index>
double Uni_B1<index>::area() const {
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

template<typename index>
double Uni_B1<index>::area(const r2degree& bound) const {
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

template<typename index>
double Uni_B1<index>::area(const pair<r2degree>& bounds) const {
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

template<typename index>
void Uni_B1<index>::compute_area_polynomial(const pair<r2degree>& bounds) {
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
    area_polynomial[3] += num_rows * range_area;

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
    for(index i = 0; i < 4; i++){
        area_polynomial[i] /= range_area;
    }
}

template<typename index>
double Uni_B1<index>::evaluate_area_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    return area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + area_polynomial[3]*x*y;
}

template<typename index>
double Uni_B1<index>::evaluate_slope_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    int k = this->d1.get_num_rows();
    return k / (area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + area_polynomial[3]*x*y);
}

template<typename index>
double Uni_B1<index>::slope() const {
    double area = this->area();
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. Consider passing a bound." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}

template<typename index>
double Uni_B1<index>::slope(const r2degree& bound) const {
    double area = this->area(bound);
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. The bound you passed is insufficient." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}

template<typename index>
double Uni_B1<index>::slope(const pair<r2degree>& bounds) const {
    double area = this->area(bounds);
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. The bound you passed is insufficient." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}


} // namespace hnf