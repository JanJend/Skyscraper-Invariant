#include "hnf_at.hpp"

namespace hnf {



bool slope_comparator::operator()(const Uni_B1& X, const Uni_B1& Y) const noexcept {
    return X.slope_value > Y.slope_value;
}

Uni_B1 find_scss_bruteforce(const R2Mat& X,
        vec<vec<vec<SparseMatrix<int>>>>& subspaces,
        R2Mat& max_subspace,
        const pair<r2degree>& bounds,
        const bool filter) {
    int k = X.get_num_rows();
    Uni_B1 scss = Uni_B1(X);
    double max_slope = scss.slope(bounds);
    if(k == 1){
        // Nothing to do?
    } else {
        assert(k < 7);
        if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
        }
        for(auto grassmanian : subspaces[k-1]){
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
            double slope = res.slope(bounds);
            if(slope > max_slope){
                max_slope = slope;
                scss = std::move(res);
                max_subspace = subspace;
            }
        }
        }
    }
    scss.slope_value = max_slope;
    return scss;
}

void skyscraper_invariant(const R2Mat& input,
    HN_factors& result,
    vec<vec<vec<SparseMatrix<int>>>>& subspaces,
    const pair<r2degree>& bounds) {
    
    R2Mat X = input;

    if(X.get_num_rows() > 1){
        
       // X.to_stream_r2(std::cout);
    }

    if(X.get_num_rows() >= 7){
        std::cout << "  Warning: Computing HNF for a module of dimension "
            << X.get_num_rows() << std::endl;
    }
    while(X.get_num_rows() > 1){
        R2Mat subspace;
        result.emplace_back(find_scss_bruteforce(X, subspaces, subspace, bounds));
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