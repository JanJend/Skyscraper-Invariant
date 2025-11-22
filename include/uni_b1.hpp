#pragma once

#ifndef UNI_B1_HPP
#define UNI_B1_HPP

#include "grlina/graded_linalg.hpp"

using namespace graded_linalg;

namespace hnf{

using R2Mat = R2GradedSparseMatrix<int>;

// Forward declarations

struct Uni_B1;
struct Slope_subdivision;

struct Uni_B1 : R2Resolution<int> {
    R2Mat d1;
    R2Mat d2;
    double slope_value;
    std::array<double, 3> area_polynomial;
    // std::unique_ptr<Slope_subdivision> slope_subdiv;  

    void compute_arrangement_quotients(vec<SparseMatrix<int>> subspaces);
    void compute_slope_subdivision(const pair<r2degree>& bounds, 
        const vec<vec<SparseMatrix<int>>>& subspaces, 
        const r2degree& cell_start, 
        const r2degree& cell_boundary);

    Uni_B1() = default;
    Uni_B1(R2Mat&& d1_, bool is_minimal = false);
    Uni_B1(const R2Mat& d1_, bool is_minimal = false);

    // Copy constructor
    Uni_B1(const Uni_B1& other) 
        : d1(other.d1), d2(other.d2), area_polynomial(other.area_polynomial) {
    }
    
    // Move constructor
    Uni_B1(Uni_B1&& other) = default;
    
    // Copy assignment
    Uni_B1& operator=(const Uni_B1& other) {
        if (this != &other) {
            d1 = other.d1;
            d2 = other.d2;
            area_polynomial = other.area_polynomial;
        }
        return *this;
    }
    
    // Move assignment
    Uni_B1& operator=(Uni_B1&& other) = default;

    double area() const;
    double area(const r2degree& bound) const;
    double area(const pair<r2degree>& bounds) const;
    void compute_area_polynomial(const pair<r2degree>& bounds);
    double evaluate_area_polynomial(r2degree d);
    double evaluate_slope_polynomial(r2degree d);
    double slope() const;
    double slope(const r2degree& bound) const;
    double slope(const pair<r2degree>& bounds) const;

    
};

} // namespace hnf

#endif // UNI_B1_HPP
