#pragma once

#ifndef HNF_AT_HEADER_HPP
#define HNF_AT_HEADER_HPP

#include "uni_b1.hpp"

namespace hnf {

using HN_factors = vec<Uni_B1>;

struct slope_comparator{
    bool operator()(const Uni_B1& X, const Uni_B1& Y) const noexcept;
};

HN_factors split_into_intervals(Uni_B1& stable_module);

HN_factors k_merge(vec<HN_factors>& sorted_lists);

void recalculate_slopes(HN_factors& composition_factors);

Uni_B1 find_scss_bruteforce(const R2Mat& X,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        R2Mat& max_subspace,
        const pair<r2degree>& bounds,
        const bool filter = false);

void skyscraper_invariant(const R2Mat& input,
    vec<HN_factors>& result,
    vec<vec<vec<SparseMatrix<int>>>>& subspaces,
    const pair<r2degree>& bounds, const bool filter = false);

template<typename Container>
vec<HN_factors> skyscraper_invariant_sum(Container& summands,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        const pair<r2degree>& bounds, const bool filter = false) {
    vec<HN_factors> result;
    for(R2Mat& X : summands){
        skyscraper_invariant(X, result, subspaces, bounds, filter);
    }
    return result;
}

template<typename Container>
void skyscraper_invariant_sum_append(Container& summands, 
        vec<HN_factors> & result,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        const pair<r2degree>& bounds, const bool filter = false) {
    for(R2Mat& X : summands){
        skyscraper_invariant(X, result, subspaces, bounds, filter);
    }
}



} // namespace hnf

#endif // HNF_AT_HEADER_HPP