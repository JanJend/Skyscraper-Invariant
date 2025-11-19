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

template<typename index>
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

template<typename index>
Point_2 dual_vertex_from_facet(Polyhedron::Facet_const_handle fit) { 
    auto he = fit->halfedge();
    Point_3 p1 = he->vertex()->point();
    Point_3 p2 = he->next()->vertex()->point();
    Point_3 p3 = he->next()->next()->vertex()->point();
    auto normal = CGAL::cross_product(p2 - p1, p3 - p1);
    
    double a = -normal.x() / normal.z();
    double b = -normal.y() / normal.z();
    return Point_2(a, b);
}

template<typename index>
bool is_upper_facet(Polyhedron::Facet_const_handle fit) {  
    auto he = fit->halfedge();
    Point_3 p1 = he->vertex()->point();
    Point_3 p2 = he->next()->vertex()->point();
    Point_3 p3 = he->next()->next()->vertex()->point();
    auto normal = CGAL::cross_product(p2 - p1, p3 - p1);
    return normal.z() > 0;
}

template<typename index>
std::optional<Segment_2> clip_segment_to_box(const Point_2& p1, const Point_2& p2, 
                                               const BoundingBox<index>& box) {
    auto result = CGAL::intersection(Segment_2(p1, p2), box.to_cgal_rect());
    if (result) {
        if (const Segment_2* seg = boost::get<Segment_2>(&*result)) {
            return *seg;
        }
    }
    return std::nullopt;
}

template<typename index>
std::optional<Segment_2> clip_ray_to_box(const Point_2& source, const K::Vector_2& dir,
                                          const BoundingBox<index>& box) {
    K::Ray_2 ray(source, dir);
    auto result = CGAL::intersection(ray, box.to_cgal_rect());
    if (result) {
        if (const Segment_2* seg = boost::get<Segment_2>(&*result)) {
            return *seg;
        }
    }
    return std::nullopt;
}

template<typename index>
std::vector<Segment_2> extract_segments_from_hull(const Polyhedron& hull, 
                                                    const BoundingBox<index>& box) {
    std::vector<Segment_2> clipped_segments;
    
    for (auto fit = hull.facets_begin(); fit != hull.facets_end(); ++fit) {
        if (!is_upper_facet<index>(fit)) continue;
        
        Point_2 arr_vertex = dual_vertex_from_facet<index>(fit);
        auto he_circ = fit->halfedge();
        
        do {
            auto opposite_facet = he_circ->opposite()->facet();
            if (opposite_facet == Polyhedron::Facet_handle()) {
                std::cerr << "Warning: Facet has no opposite facet." << std::endl;
                he_circ = he_circ->next();
                continue;
            }
            
            if (is_upper_facet<index>(opposite_facet)) {
                // Edge between two lower facets
                Point_2 arr_vertex_opp = dual_vertex_from_facet<index>(opposite_facet);
                if (auto seg = clip_segment_to_box<index>(arr_vertex, arr_vertex_opp, box)) {
                    clipped_segments.push_back(*seg);
                }
            } else {
                // Ray case
                Point_3 edge_start = he_circ->vertex()->point();
                Point_3 edge_end = he_circ->next()->vertex()->point();
                K::Vector_2 edge_dir(edge_end.x() - edge_start.x(),
                                    edge_end.y() - edge_start.y());
                K::Vector_2 ray_dir(-edge_dir.y(), edge_dir.x());
                
                auto he_opp = opposite_facet->halfedge();
                Point_3 op1 = he_opp->vertex()->point();
                Point_3 op2 = he_opp->next()->vertex()->point();
                Point_3 op3 = he_opp->next()->next()->vertex()->point();
                auto normal_opp = CGAL::cross_product(op2 - op1, op3 - op1);
                
                if (normal_opp.z() < 0) {
                    ray_dir = -ray_dir;
                }
                
                if (auto seg = clip_ray_to_box<index>(arr_vertex, ray_dir, box)) {
                    clipped_segments.push_back(*seg);
                }
            }
            he_circ = he_circ->next();
        } while (he_circ != fit->halfedge());
    }
    
    return clipped_segments;
}

template<typename index>
void assign_face_data(Arrangement<index>& arr, const Polyhedron& hull,
                     const std::vector<std::array<double,3>>& polynomials,
                     const std::map<Point_3, size_t>& point_to_index) {
    CGAL::Arr_naive_point_location<Arrangement<index>> pl(arr);
    
    for (auto vit = hull.vertices_begin(); vit != hull.vertices_end(); ++vit) {
        Point_3 hull_vertex = vit->point();
        
        // Check if this vertex is on the lower hull
        bool is_lower = false;
        auto he_circ = vit->halfedge();
        do {
            if (he_circ->facet() != Polyhedron::Facet_handle() && 
                is_upper_facet<index>(he_circ->facet())) {
                is_lower = true;
                break;
            }
            he_circ = he_circ->next()->opposite();
        } while (he_circ != vit->halfedge());
        
        if (!is_lower) continue;
        
        Point_2 query(hull_vertex.x(), hull_vertex.y());
        auto obj = pl.locate(query);
        
        typename Arrangement<index>::Face_handle fh;
        if (CGAL::assign(fh, obj) && !fh->is_unbounded()) {
            size_t orig_idx = point_to_index.at(hull_vertex);
            fh->data().subspace_index = orig_idx;
            fh->data().slope_polynomial = polynomials[orig_idx];
        }
    }
}

template<typename index>
Arrangement<index> subdivision_from_polynomials(const vec<std::array<double,3>>& polynomials,
                                                 const r2degree& cell_start, 
                                                 const r2degree& cell_end) {
    BoundingBox<index> box{cell_start.first, cell_end.first, 
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
    std::vector<Segment_2> clipped_segments = extract_segments_from_hull<index>(hull, box);
    
    // Build arrangement with bounding box
    Arrangement<index> arr;
    auto bbox_segments = box.boundary_segments();
    CGAL::insert(arr, bbox_segments.begin(), bbox_segments.end());
    CGAL::insert(arr, clipped_segments.begin(), clipped_segments.end());
    
    // Assign face data
    assign_face_data(arr, hull, polynomials, point_to_index);
    
    return arr;
}

template<typename index>
void Uni_B1<index>::compute_arrangement_quotients(vec<SparseMatrix<index>> subspaces){
    //TO-DO implement
};

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
    vec<std::array<double, 3>> slope_polynomials;
    if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
    }

    for(size_t i = 1; i < subspaces[k-1].size(); i++){
        //skip i = 0, because it is the empty space
        auto ungraded_subspace = subspaces[k-1][i];
        int num_gens = ungraded_subspace.get_num_cols();
        R2Mat subspace = R2Mat(ungraded_subspace);
        subspace.row_degrees = X.row_degrees;
        subspace.col_degrees = vec<r2degree>(num_gens, X.row_degrees[0]);
            assert(subspace.get_num_rows() == X.get_num_rows());
            assert(subspace.get_num_cols() == num_gens);
        R2Mat submodule_pres = X.submodule_generated_by(subspace);
        X.column_reduction_graded(); //full minimisation should not be necessary
        Uni_B1<int> res(submodule_pres);
        res.compute_area_polynomial(bounds);  // Compute first
        slope_polynomials.emplace_back(res.area_polynomial);
        for(auto& coeff : slope_polynomials.back()){
            coeff /= static_cast<double>(num_gens);
        }
    }

    this->slope_subdiv = std::make_unique<Slope_subdivision<index>>(
        subdivision_from_polynomials<index>(slope_polynomials, X.row_degrees[0], cell_boundary)
    );

    this->compute_arrangement_quotients(subspaces[k-1]);
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
    for(index i = 0; i < 3; i++){
        area_polynomial[i] /= range_area;
    }
}

template<typename index>
double Uni_B1<index>::evaluate_area_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    return area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + x*y;
}

template<typename index>
double Uni_B1<index>::evaluate_slope_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    int k = this->d1.get_num_rows();
    return k / (area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + x*y);
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

// Explicit template instantiation for int
template class hnf::Uni_B1<int>;
template hnf::Arrangement<int> hnf::subdivision_from_polynomials(
    const vec<std::array<double,3>>&, const r2degree&, const r2degree&);

} // namespace hnf