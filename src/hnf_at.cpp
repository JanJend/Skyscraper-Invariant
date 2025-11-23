#include "hnf_at.hpp"

namespace hnf {


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
        assert(k < 7);
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
    HN_factors& result,
    vec<vec<vec<SparseMatrix<int>>>>& subspaces,
    const pair<r2degree>& bounds,
    const bool filter) {
    
    R2Mat X = input;

    if(X.get_num_rows() > 1){
        
       // X.to_stream_r2(std::cout);
    }

    if(X.get_num_rows() >= 7){
        std::cout << "  Warning: Computing HNF for a module of dimension "
            << X.get_num_rows() << std::endl;
    }
    result.reserve(X.get_num_rows());
    while(X.get_num_rows() > 1){
        R2Mat subspace;
        result.emplace_back(find_scss_bruteforce(X, subspaces, subspace, bounds, filter));
        if(result.back().d1.get_num_rows() == X.get_num_rows()){
            break;
        } else {
            X.quotient_by(subspace);
        }
    }
    if(X.get_num_rows() == 1){
        Uni_B1 res(X);
        res.slope_value = res.slope(bounds);
        result.emplace_back(std::move(res));
    } else if (X.get_num_rows() == 0){
        // Nothing to do.
    }
}


} // namespace hnf