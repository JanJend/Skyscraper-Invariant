#pragma once

#ifndef UNI_B1_HPP
#define UNI_B1_HPP

#include "dynamic_matrix_adapter.hpp"

using namespace graded_linalg;

namespace hnf{

inline R2Mat& module_presentation(R2Mat& matrix) { return matrix; }
inline const R2Mat& module_presentation(const R2Mat& matrix) { return matrix; }
// Invariant algorithms reorder columns and compute scratch grid state. Give
// them a working copy so querying a module preserves every resolution map.
inline R2Mat module_presentation(const R2PModule& module) {
    R2Mat working = module.presentation();
    return working;
}
inline R2Mat module_presentation(R2PModule& module) {
    return module_presentation(static_cast<const R2PModule&>(module));
}
template<class Input, std::enable_if_t<needs_editable_presentation_v<Input>, int> = 0>
R2Mat module_presentation(const Input& input) { return editable_presentation(input); }
template<class Input, std::enable_if_t<needs_editable_presentation_v<Input>, int> = 0>
R2Mat module_presentation(Input& input) { return editable_presentation(input); }

template <typename MatrixLike, std::enable_if_t<std::is_base_of_v<R2Mat, MatrixLike>, int> = 0>
R2Mat& module_presentation(MatrixLike& matrix) { return static_cast<R2Mat&>(matrix); }

template <typename MatrixLike, std::enable_if_t<std::is_base_of_v<R2Mat, MatrixLike>, int> = 0>
const R2Mat& module_presentation(const MatrixLike& matrix) { return static_cast<const R2Mat&>(matrix); }

// Forward declarations

struct Uni_B1;

struct Uni_B1 {
    R2Mat d1;
    R2Mat d2;
    double slope_value;
    std::array<double, 3> area_polynomial;

    void compute_arrangement_quotients(vec<SparseMatrix<int>> subspaces);
    void compute_slope_subdivision(const pair<r2degree>& bounds, 
        const vec<vec<SparseMatrix<int>>>& subspaces, 
        const r2degree& cell_start, 
        const r2degree& cell_boundary);

    Uni_B1() = default;
    Uni_B1(R2Mat&& d1_, bool is_minimal = false);
    Uni_B1(const R2Mat& d1_, bool is_minimal = false);
    Uni_B1(const R2PModule& module, bool is_minimal = false);
    template<class Matrix, std::enable_if_t<is_adapted_matrix<Matrix>::value, int> = 0>
    Uni_B1(const Matrix& matrix, bool is_minimal = false)
        : Uni_B1(editable_presentation(matrix), is_minimal) {}
    template<class Matrix, std::enable_if_t<is_adapted_matrix<Matrix>::value, int> = 0>
    Uni_B1(const Module<Matrix>& module, bool is_minimal = false)
        : d1(editable_presentation(module)) {
        if (module.proj_res_.size() > 1)
            d2 = editable_presentation(module.proj_res_[1]);
        else
            d2 = Uni_B1(d1, is_minimal).d2;
    }

    // Copy constructor
    Uni_B1(const Uni_B1& other) 
        : d1(other.d1), d2(other.d2), area_polynomial(other.area_polynomial), slope_value(other.slope_value) {
    }
    
    // Move constructor (already have)
    Uni_B1(Uni_B1&& other) = default;

    // Move assignment
    Uni_B1& operator=(Uni_B1&& other) = default;
    
    // Copy assignment
    Uni_B1& operator=(const Uni_B1& other) {
        if (this != &other) {
            d1 = other.d1;
            d2 = other.d2;
            area_polynomial = other.area_polynomial;
            slope_value = other.slope_value;
        }
        return *this;
    }


    double area() const;
    double area(const r2degree& bound) const;
    double area_unbounded(const pair<r2degree>& bounds) const;
    double area(const pair<r2degree>& bounds) const;
    void compute_area_polynomial(const pair<r2degree>& bounds);
    void compute_area_polynomial_unbounded(const pair<r2degree>& bounds);
    double evaluate_area_polynomial(r2degree d, const pair<r2degree>& bounds);
    double evaluate_slope_polynomial(r2degree d, const pair<r2degree>& bounds);
    double slope() const;
    double slope(const r2degree& bound) const;
    double slope(const pair<r2degree>& bounds) const;
    void compute_slope(const pair<r2degree>& bounds);
    
};

} // namespace hnf

#endif // UNI_B1_HPP
