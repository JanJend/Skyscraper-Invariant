#include "hnf_at.hpp"
#include <iostream>
#include <filesystem>

using namespace graded_linalg;
using namespace hnf;

void hnf_at_origin(std::filesystem::path input_path) {
    
    R2GradedSparseMatrix<int> X = R2GradedSparseMatrix<int>(input_path.string());
    int dim = X.get_num_rows();
    pair<r2degree> bounds = X.bounding_box();
    vec<vec<vec<SparseMatrix<int>>>> subspaces = sparse_seperated_grassmannians<int>(dim);
    HN_factors result;
    skyscraper_invariant(X, result, subspaces, bounds);
    for(auto& factor : result){
        std::cout << "Slope: " << factor.slope_value << " Module: \n";
        factor.d1.to_stream_r2(std::cout);
    }
}


int main(int argc, char** argv) {
    
    std::string filepath;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file_path>" << std::endl;
        return 1;
    } else {
        filepath = argv[1];
    }

    std::filesystem::path input_path(filepath);
    
    hnf_at_origin(input_path);
    return 0;
} // main