#ifndef HNF_HEADER_HPP
#define HNF_HEADER_HPP

#include "aida_interface.hpp"
#include "grlina/graded_linalg.hpp"
#include <iostream>
#include <unistd.h>
#include <getopt.h>

using namespace graded_linalg;

namespace hnf{

template<typename index>
struct Uni_B1{
    R2GradedSparseMatrix<index> d1;
    R2GradedSparseMatrix<index> d2;
    std::array<double, 4> area_polynomial;

    Uni_B1() = default;
    Uni_B1(R2GradedSparseMatrix<index>&& d1_, bool is_minimal = false);
    Uni_B1(const R2GradedSparseMatrix<index>& d1_, bool is_minimal = false);

    index dim_at(r2degree alpha) const;
    double area() const;
    double area(const r2degree& bound) const;
    double area(const pair<r2degree>& bounds) const;
    void compute_area_polynomial(const pair<r2degree>& bounds);
    double evaluate_area_polynomial(r2degree d);
    double evaluate_slope_polynomial(r2degree d);
    double slope() const;
    double slope(const r2degree& bound) const;
    double slope(const pair<r2degree>& bounds) const;
    array<index> dimension_vector(bool sort = true) const;
};

using Block = aida::Block;
using Module_w_slope = std::pair<Uni_B1<int>, double>;
using Block_list = aida::Block_list;
using HN_factors = vec<Module_w_slope>;
using R2Mat = R2GradedSparseMatrix<int>;

struct slope_comparator{
    bool operator()(const Module_w_slope& X, const Module_w_slope& Y) const noexcept;
};

Module_w_slope find_scss_bruteforce(const R2Mat& X,
        vec<vec<SparseMatrix<int>>>& subspaces,
        R2Mat& max_subspace,
        const pair<r2degree>& bounds);

HN_factors skyscraper_invariant(Block_list& summands,
        vec<vec<SparseMatrix<int>>>& subspaces,
        const pair<r2degree>& bounds);

void skyscraper_invariant(const R2Mat& input,
    HN_factors& result,
    vec<vec<SparseMatrix<int>>>& subspaces,
    const pair<r2degree>& bounds);

void calculate_stats(const std::vector<int>& all_dimensions);

vec<r2degree> get_grid_diagonal(pair<r2degree> bounds, int grid_length);

r2degree get_grid_step(const r2degree& lower_bound, const r2degree& upper_bound,
    const int& grid_length_x, const int& grid_length_y);

template<typename Outputstream>
void to_stream(Outputstream& ostream, Module_w_slope& scss);

void write_slopes_to_csv(const vec<vec<double>>& slopes,
        const vec<r2degree>& grid_points,
        const std::string& filename);

void show_progress_bar(int& i, int& total, std::string& name);

template<typename Container, typename Outputstream>
void process_summands_fixed_grid(aida::AIDA_functor& decomposer,
    Outputstream& ostream,
    const int& grid_length_x, const int& grid_length_y,
    Container& indecomps);

struct Dynamic_HNF {
    vec<vec<Uni_B1<int>>> indecomposable_summands;
    vec<int> grid_ind_dimensions;

    Dynamic_HNF();
    void compute_HNF_row(aida::AIDA_functor& decomposer,
        R2Mat& M,
        int& y_index,
        pair<r2degree> slope_bounds);
};

template<typename Container>
std::tuple<r2degree, r2degree, r2degree, pair<r2degree>> compute_bounds_and_grid(
    Container& indecomps,
    vec<int>& first_ind_dimensions,
    const int& grid_length_x,
    const int& grid_length_y);

template<typename Container>
void update_HNF_rows_at_y_level(
    r2degree current_grid_degree,
    Container& indecomps,
    vec<pair<int>>& grid_locations,
    vec<Dynamic_HNF>& local_grid_row_data,
    aida::AIDA_functor& decomposer,
    const pair<r2degree>& slope_bounds);

template<typename Container>
void update_grid_locations_x(
    r2degree current_grid_degree,
    Container& indecomps,
    vec<pair<int>>& grid_locations);

template<typename Outputstream>
void write_grid_metadata(Outputstream& ostream,
    int grid_length_x, int grid_length_y,
    const r2degree& lower_bound,
    const r2degree& upper_bound,
    const r2degree& grid_step,
    const pair<r2degree>& slope_bounds,
    bool show_info = false);

bool essentially_equal(double a, double b, double relTol = 1e-9, double absTol = 1e-12);

void compare_slopes_test(
    r2degree current_grid_degree,
    r2degree local_grid_degree,
    const HN_factors& composition_factors,
    const HN_factors& test_factors,
    int i, int j, int k);

template<typename Container, typename Outputstream>
void process_grid_cell(
    int i, int j,
    r2degree current_grid_degree,
    Container& indecomps,
    vec<pair<int>>& grid_locations,
    vec<Dynamic_HNF>& local_grid_row_data,
    array<HN_factors>& composition_factors,
    vec<int>& grid_ind_dimensions,
    vec<int>& all_scss_dimensions,
    Outputstream& ostream,
    vec<vec<SparseMatrix<int>>>& subspaces,
    const pair<r2degree>& slope_bounds,
    aida::AIDA_functor& decomposer);

template<typename Container, typename Outputstream>
void process_summands_smart_grid(aida::AIDA_functor& decomposer,
    Outputstream& ostream,
    const int& grid_length_x, const int& grid_length_y,
    Container& indecomps);

template<typename Outputstream>
void full_grid_induced_decomposition(aida::AIDA_functor& decomposer,
    std::ifstream& istream, Outputstream& ostream,
    bool show_indecomp_statistics, bool show_runtime_statistics,
    bool dynamic_grid = true,
    bool is_decomposed = false,
    const int& grid_length_x = 30, const int& grid_length_y = 30);

} // namespace hnf


#endif