#include "hnf_at.hpp"
#include <iostream>
#include <filesystem>

using namespace graded_linalg;
using namespace hnf;

void hnf_at_origin(std::filesystem::path input_path, r2degree upper_bound, bool filter) {
    
    R2GradedSparseMatrix<int> X = R2GradedSparseMatrix<int>(input_path.string());
    int dim = X.get_num_rows();
    pair<r2degree> bounds = X.bounding_box();
    if(upper_bound != r2degree{0,0}){
        bounds.second = upper_bound;
    } else {
        double padding = 0.1;
        r2degree range = bounds.second - bounds.first;
        bounds.second = bounds.second + r2degree{padding*range.first, padding*range.second};
    }
    vec<vec<vec<SparseMatrix<int>>>> subspaces = sparse_seperated_grassmannians<int>(dim);
    HN_factors result;
    skyscraper_invariant(X, result, subspaces, bounds, filter);
    for(auto& factor : result){
        std::cout << "Slope: " << factor.slope_value << " Module: \n";
        factor.d1.to_stream_r2(std::cout);
    }
}


int main(int argc, char** argv) {
    
    std::string filepath;
    r2degree upper_bound = {0,0};
    bool filter = false;
    if (argc < 2 || argc > 4) {
        std::cerr << "Usage: " << argv[0] << " <file_path> <upper_bound> <filter>" << std::endl;
        return 1;
        //    filepath = "/home/wsljan/MP-Workspace/data/hypoxic_regions/hypoxic2_FoxP3_dim1_200x200_snapped_induced/4_0.scc";
        //    upper_bound = {0.5, 2.2};
        //    filter = true;
    } else {
        filepath = argv[1];
        if (argc >= 3) {
            std::istringstream bound_stream(argv[2]);
            char ignore;
            double x, y;
            bound_stream >> ignore >> x >> ignore >> y >> ignore;
            upper_bound = {x, y};
        }
        if (argc == 4) {
            std::istringstream filter_stream(argv[3]);
            filter_stream >> std::boolalpha >> filter;
        }
    }

    std::filesystem::path input_path(filepath);
    
    hnf_at_origin(input_path, upper_bound, filter);
    return 0;
} // main