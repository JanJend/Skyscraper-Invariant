#include "uni_b1.hpp"

namespace hnf {

Uni_B1::Uni_B1(R2Mat&& d1_, bool is_minimal)
    : R2Resolution<int>(std::move(d1_)) {
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
    : R2Resolution<int>(d1_) {
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