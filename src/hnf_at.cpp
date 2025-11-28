#include "hnf_at.hpp"

namespace hnf {

HN_factors split_into_intervals(Uni_B1& stable_module){
    R2Mat pres = stable_module.d1;
    HN_factors intervals;
    double slope = stable_module.slope_value;
    
    while(pres.get_num_rows() > 0){
        R2Mat subspace = R2Mat(1, pres.get_num_rows());
        subspace.data = {{0}};
        subspace.row_degrees = pres.row_degrees;
        subspace.col_degrees = vec<r2degree>(1, pres.row_degrees[0]);
        R2Mat submodule_pres = pres.submodule_generated_by(subspace);
        submodule_pres.minimize();
        intervals.emplace_back(Uni_B1(submodule_pres));
        intervals.back().slope_value = slope;
        pres.quotient_by(subspace);
    }
    
    return intervals;
}

vec<Uni_B1>  k_merge(vec<vec<Uni_B1>>& factors_from_indecomp){
    vec<Uni_B1> result;
    // Max-heap: pair of (slope_value, vector_index)
    auto comp = [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
        return a.first < b.first;
    };
    std::priority_queue<std::pair<double, int>, 
                       vec<std::pair<double, int>>, 
                       decltype(comp)> pq(comp);
    
    // Track current position in each vector
    vec<int> indices(factors_from_indecomp.size(), 0);
    
    // Initialize heap with slope_value of first element in each vector
    for (size_t i = 0; i < factors_from_indecomp.size(); ++i) {
        if (!factors_from_indecomp[i].empty()) {
            pq.push({factors_from_indecomp[i][0].slope_value, static_cast<int>(i)});
        }
    }
    
    // Extract max, move object, advance that vector
    while (!pq.empty()) {
        auto [slope_value, vec_idx] = pq.top();
        pq.pop();
        
        result.push_back(std::move(factors_from_indecomp[vec_idx][indices[vec_idx]]));
        indices[vec_idx]++;
        
        if (indices[vec_idx] < factors_from_indecomp[vec_idx].size()) {
            pq.push({factors_from_indecomp[vec_idx][indices[vec_idx]].slope_value, vec_idx});
        }
    }
    
    return result;
}


void recalculate_slopes(HN_factors& composition_factors) {
    if(composition_factors.empty()){
        return;
    }
    double last_slope = composition_factors.front().slope_value;
    int aggregated_dimension = composition_factors.front().d1.get_num_rows();
    for(size_t i = 1; i < composition_factors.size(); i++){
        double current_slope = composition_factors[i].slope_value;
        int current_dimension = composition_factors[i].d1.get_num_rows();
        double area_last = static_cast<double>(aggregated_dimension) / last_slope; 
        double area_current = static_cast<double>(current_dimension) / current_slope;
        double total_area = area_last + area_current;
        double new_slope = static_cast<double>(aggregated_dimension + current_dimension) / total_area;
        composition_factors[i].slope_value = new_slope;
        last_slope = new_slope;
        aggregated_dimension += current_dimension;
    }
}

bool slope_comparator::operator()(const Uni_B1& X, const Uni_B1& Y) const noexcept {
    return X.slope_value > Y.slope_value;
}
vec<Uni_B1> get_slopes_at_dim(const R2Mat& X,
    vec<SparseMatrix<int>>& grassmanian,
    Uni_B1& scss,
    R2Mat& max_subspace,
    const pair<r2degree>& bounds) {
    vec<Uni_B1> submodules;
    for(auto ungraded_subspace : grassmanian){
        int num_gens = ungraded_subspace.get_num_cols();
        if(num_gens == 0){
            continue;
        }
        R2Mat subspace = R2Mat(ungraded_subspace);
        subspace.row_degrees = X.row_degrees;
        subspace.col_degrees = vec<r2degree>(num_gens, X.row_degrees[0]);
        assert(subspace.get_num_rows() == X.get_num_rows());
        assert(subspace.get_num_cols() == num_gens);
        R2Mat submodule_pres = X.submodule_generated_by(subspace);
        submodule_pres.minimize();
        submodules.emplace_back(Uni_B1(submodule_pres));
        submodules.back().slope_value = submodules.back().slope(bounds);
        if(submodules.back().slope_value > scss.slope_value){
            scss = submodules.back();
            max_subspace = subspace;
        }
    }
    return submodules;
}

void find_scss_of_dim(const R2Mat& X,
    vec<SparseMatrix<int>>& grassmanian,
    Uni_B1& scss,
    R2Mat& max_subspace,
    const pair<r2degree>& bounds) {
    for(auto ungraded_subspace : grassmanian){
        int num_gens = ungraded_subspace.get_num_cols();
        if(num_gens == 0){
            continue;
        }
        R2Mat subspace = R2Mat(ungraded_subspace);
        subspace.row_degrees = X.row_degrees;
        subspace.col_degrees = vec<r2degree>(num_gens, X.row_degrees[0]);
        assert(subspace.get_num_rows() == X.get_num_rows());
        assert(subspace.get_num_cols() == num_gens);
        R2Mat submodule_pres = X.submodule_generated_by(subspace);
        submodule_pres.minimize();
        Uni_B1 res(submodule_pres);
        res.slope_value = res.slope(bounds);
        if(res.slope_value > scss.slope_value){
            scss = res;
            max_subspace = subspace;
        }
    }
}

Uni_B1 find_scss_bruteforce(const R2Mat& X,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        R2Mat& max_subspace,
        const pair<r2degree>& bounds,
        const bool filter) {
    int k = X.get_num_rows();
    Uni_B1 scss = Uni_B1(X);
    scss.slope_value = 0.0;
    if(k == 1){
        // Nothing to do?
    } else {
        assert(k < 8);
        if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
        }
        vec<bool> filtered_out(k+1, false);
        if(filter){
            vec<Uni_B1> dim_1_submodules = get_slopes_at_dim(X, subspaces[k-1][1], scss, max_subspace, bounds);
            filtered_out[1] = true;
            for(size_t j = 2; j < k+1; j++){
                size_t high_slope_count = 0;
                for(const Uni_B1& submodule : dim_1_submodules){
                    if(submodule.slope_value > scss.slope_value / static_cast<double>(j)){
                        high_slope_count++;
                    } 
                }
                if(high_slope_count < std::pow(2, j) - 1){
                    filtered_out[j] = true;
                }
            }
        }
        for(size_t i = 1; i < subspaces[k-1].size(); i++){
            if(filter){
                if(filtered_out[i]){
                    continue;
                }
            }
            auto& grassmanian = subspaces[k-1][i];
            find_scss_of_dim(X, grassmanian, scss, max_subspace, bounds);
        }
        if (filter){
            std::cout << " Filtered out all subspaces of dimensions: ";
            for(size_t i = 2; i < filtered_out.size(); i++){
                if(filtered_out[i]){
                    std::cout << i << " ";
                }
            }
            std::cout << std::endl;
        }
    }
    return scss;
}

void skyscraper_invariant(const R2Mat& input,
    vec<HN_factors>& result,
    vec<vec<vec<SparseMatrix<int>>>>& subspaces,
    const pair<r2degree>& bounds,
    const bool filter) {
    
    R2Mat X = input;

    assert(X.get_num_rows() > 1);

    if(X.get_num_rows() > 1){
        
       // X.to_stream_r2(std::cout);
    }
    if(X.get_num_rows() >= 7){
        std::cout << "  Warning: Computing HNF for a module of dimension "
            << X.get_num_rows() << std::endl;
    }

    result.emplace_back( HN_factors() );
    result.reserve(X.get_num_rows());
    while(X.get_num_rows() > 0){
        R2Mat subspace;
        result.back().emplace_back(find_scss_bruteforce(X, subspaces, subspace, bounds, filter));
        if(result.back().back().d1.get_num_rows() == X.get_num_rows()){
            break;
        } else {
            X.quotient_by(subspace);
        }
    }
}


} // namespace hnf