#include "uni_b1.hpp"

namespace hnf {

Uni_B1::Uni_B1(R2Mat&& d1_, bool is_minimal)
    : Uni_B1(static_cast<const R2Mat&>(d1_), is_minimal) {}

Uni_B1::Uni_B1(const R2Mat& d1_, bool /*is_minimal*/) : d1(d1_) {
    require_two_parameters(d1);
    auto kernel_source = d1;
    d2 = kernel_source.graded_kernel();
}

Uni_B1::Uni_B1(const R2PModule& module, bool is_minimal) {
    require_two_parameters(module.presentation());
    if (module.proj_res_.size() > 1) {
        d1 = module.presentation();
        d2 = module.proj_res_[1];
    } else {
        const Uni_B1 resolution(module.presentation(), is_minimal);
        d1 = resolution.d1;
        d2 = resolution.d2;
    }
}


double Uni_B1::area() const {
    auto [min, max] = geometry_bounds(d1);
    double base_area = d1.get_num_rows()*(max[0] - min[0]) * (max[1] - min[1]);
    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        base_area -= (max[0] - degree[0]) * (max[1] - degree[1]);
    }
    for(int column = 0; column < d2.get_num_cols(); ++column){
        const auto degree = real_col_degree(d2, column);
        base_area += (max[0] - degree[0]) * (max[1] - degree[1]);
    }
    return base_area;
}


double Uni_B1::area(const r2degree& bound) const {
    auto [min, max] = geometry_bounds(d1);
    max = bound;
    double base_area = d1.get_num_rows()*(max[0] - min[0]) * (max[1] - min[1]);
    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        base_area -= (max[0] - degree[0]) * (max[1] - degree[1]);
    }
    for(int column = 0; column < d2.get_num_cols(); ++column){
        const auto degree = real_col_degree(d2, column);
        base_area += (max[0] - degree[0]) * (max[1] - degree[1]);
    }
    return base_area;
}

double Uni_B1::area(const pair<r2degree>& bounds) const {
    r2degree range = bounds.second-bounds.first;
    double range_area = range[0] * range[1];
    double area = 0;
    area += d1.get_num_rows()*real_row_degree(d1, 0)[0] * real_row_degree(d1, 0)[1];
    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        area -= degree[0] * degree[1];
    }
    for(int column = 0; column < d2.get_num_cols(); ++column){
        const auto degree = real_col_degree(d2, column);
        area += degree[0] * degree[1];
    }
    double normalised_area = area/range_area;
    return normalised_area;
}


double Uni_B1::area_unbounded(const pair<r2degree>& bounds) const {
    auto [min, max] = geometry_bounds(d1);
    max = bounds.second;
    assert(max[0] >= min[0]);
    assert(max[1] >= min[1]);
    assert(min == real_row_degree(d1, 0));
    r2degree range = bounds.second-bounds.first;
    double range_area = range[0] * range[1];
    double base_area = d1.get_num_rows()*(max[0] - min[0])*(max[1] - min[1]);
    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        assert(degree[0] <= max[0]);
        assert(degree[1] <= max[1]);
        assert(degree[0] >= min[0]);
        assert(degree[1] >= min[1]);
        base_area -= (max[0] - degree[0]) * (max[1] - degree[1]);
    }
    for(int column = 0; column < d2.get_num_cols(); ++column){
        const auto degree = real_col_degree(d2, column);
        assert(degree[0] <= max[0]);
        assert(degree[1] <= max[1]);
        assert(degree[0] >= min[0]);
        assert(degree[1] >= min[1]);
        base_area += (max[0] - degree[0]) * (max[1] - degree[1]);
    }
    double normalised_area = base_area/range_area;
    return normalised_area;
}

void Uni_B1::compute_area_polynomial(const pair<r2degree>& bounds) {

    int num_rows = d1.get_num_rows();
    r2degree range = bounds.second-bounds.first;
    double range_area = range[0] * range[1];

    area_polynomial.fill(0.0);
    area_polynomial[0] = this->area(bounds);
    area_polynomial[1] += num_rows * real_row_degree(d1, 0)[1];
    area_polynomial[2] += num_rows * real_row_degree(d1, 0)[0];
    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        if(degree[0] == real_row_degree(d1, 0)[0]){
            area_polynomial[1] -= degree[1];
        }
        if(degree[1] == real_row_degree(d1, 0)[1]){
            area_polynomial[2] -= degree[0];
        }
    }
    for(int i = 1; i < 3; i++){
        area_polynomial[i] /= range_area;
    }

}

void Uni_B1::compute_area_polynomial_unbounded(const pair<r2degree>& bounds) {
    int num_rows = d1.get_num_rows();
    assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(),
                         Degree_traits<R2Mat::degree_type>::lex_lambda()));
    bool test = true;

    r2degree range = bounds.second-bounds.first;
    double range_area = range[0] * range[1];

    const r2degree gen_degree = real_row_degree(d1, 0);
    area_polynomial.fill(0.0);

    r2degree gen_vector = bounds.second - gen_degree;
    area_polynomial[0] += num_rows * gen_vector[0] * gen_vector[1];
    area_polynomial[1] -= num_rows * gen_vector[1];
    area_polynomial[2] -= num_rows * gen_vector[0];

    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        assert( Degree_traits<r2degree>::smaller_equal(degree, bounds.second));
        assert( Degree_traits<r2degree>::smaller_equal(bounds.first, degree) );
        r2degree rel_vector = bounds.second - degree;
        area_polynomial[0] -= rel_vector[0]* rel_vector[1];
    }

    for(int column = 0; column < d1.get_num_cols(); ++column){
        const auto degree = real_col_degree(d1, column);
        if(degree[0] == gen_degree[0]){
            assert(bounds.second[1] - degree[1] >= 0);
            area_polynomial[1] += (bounds.second[1] - degree[1]);
            if( area_polynomial[1] > 0){
                std::cout << "Area polynomial[1] is positive: " << area_polynomial[1] << std::endl;
                this->d1.print_graded();
                std::cout << "Bounds: " << bounds.first[0] << " " << bounds.first[1] << " to "
                    << bounds.second[0] << " " << bounds.second[1] << std::endl;
                std::cout << "Range area: " << range_area << std::endl;
                assert(false);
            }
        }
        if(degree[1] == gen_degree[1]){
            assert(bounds.second[0] - degree[0] >= 0);
            area_polynomial[2] += (bounds.second[0] - degree[0]);
            if( test && area_polynomial[2] > 0){
                std::cout << "Area polynomial[2] is positive: " << area_polynomial[1] << std::endl;
                this->d1.print_graded();
                std::cout << "Bounds: " << bounds.first[0] << " " << bounds.first[1] << " to "
                    << bounds.second[0] << " " << bounds.second[1] << std::endl;
                std::cout << "Range area: " << range_area << std::endl;
                assert(false);
            }
        }
    }

    for(int column = 0; column < d2.get_num_cols(); ++column){
        const auto degree = real_col_degree(d2, column);
        assert( Degree_traits<r2degree>::smaller_equal(degree, bounds.second));
        assert( Degree_traits<r2degree>::smaller_equal(bounds.first, degree) );
        r2degree rel_vector = bounds.second - degree;
        area_polynomial[0] += rel_vector[0]* rel_vector[1];
    }
    for(int i = 0; i < 3; i++){
        area_polynomial[i] /= range_area;
    }
}

double Uni_B1::evaluate_area_polynomial(r2degree d, const pair<r2degree>& bounds) {
    double range_area = (bounds.second[0] - bounds.first[0]) * (bounds.second[1] - bounds.first[1]);
    double& x = d[0];
    double& y = d[1];
    return area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + x*y/range_area;
}

double Uni_B1::evaluate_slope_polynomial(r2degree d, const pair<r2degree>& bounds) {
    double range_area = (bounds.second[0] - bounds.first[0]) * (bounds.second[1] - bounds.first[1]);
    double& x = d[0];
    double& y = d[1];
    double k = this->d1.get_num_rows();
    return k / (area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + x*y/range_area);
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

void Uni_B1::compute_slope(const pair<r2degree>& bounds){
    slope_value = this->slope(bounds);
}

} // namespace hnf
