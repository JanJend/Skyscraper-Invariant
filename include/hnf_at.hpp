#pragma once

#ifndef HNF_AT_HEADER_HPP
#define HNF_AT_HEADER_HPP

#include "uni_b1.hpp"

namespace hnf {

using HN_factors = vec<Uni_B1>;

struct slope_comparator{
    bool operator()(const Uni_B1& X, const Uni_B1& Y) const noexcept;
};

Uni_B1 find_scss_bruteforce(const R2Mat& X,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        R2Mat& max_subspace,
        const pair<r2degree>& bounds,
        const bool filter = false);

void skyscraper_invariant(const R2Mat& input,
    HN_factors& result,
    vec<vec<vec<SparseMatrix<int>>>>& subspaces,
    const pair<r2degree>& bounds, const bool filter = false);

template<typename Container>
HN_factors skyscraper_invariant_sum(Container& summands,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        const pair<r2degree>& bounds, const bool filter = false) {
    HN_factors result;
    for(R2Mat& X : summands){
        skyscraper_invariant(X, result, subspaces, bounds, filter);
    }
    return result;
}



} // namespace hnf

#endif // HNF_AT_HEADER_HPP