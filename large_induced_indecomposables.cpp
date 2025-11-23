#include "grlina/r2graded_matrix.hpp"
#include "grlina/graded_linalg.hpp"
#include "aida_interface.hpp"
#include <iostream>
#include <filesystem>
  
using namespace graded_linalg;




void large_induced_indecomp_submodules(R2GradedSparseMatrix<int>& pres, std::filesystem::path output_dir, vec<int> counter = vec<int>(50,0)) {
    aida::AIDA_functor decomposer = aida::AIDA_functor();
    pres.compute_col_batches();
    pres.compute_grid_representation();
    bool interval_enabled = false;
    int interval = 10;
    int pause_counter_x = interval;
    int pause_counter_y = interval;
    bool pause = false;
    vec<pair<int>> checked_list = vec<pair<int>>();
    for(int i = 0; i < pres.x_grid.size(); i++){
        const double& d = pres.x_grid[i];
        for(int j = 0; j < pres.y_grid.size(); j++){
            const double& e = pres.y_grid[j];
            std::pair<int, int> current_index = {i, j};
            bool in_interval = false;
            if(interval_enabled){
                for(const auto& [a, b] : checked_list) {
                    if(i >= a && i <= a + interval && 
                    j >= b && j <= b + interval) {
                        in_interval = true;
                        break;
                    }
                }
                if(in_interval) {
                    continue;
                }
            }
            r2degree degree = {d, e};
            R2GradedSparseMatrix<int> submodule = pres.submodule_generated_at(degree);
            submodule.compute_col_batches();
            aida::Block_list B_list;
            decomposer(submodule, B_list);
            for(const auto& ind : B_list){
                int dim = ind.get_num_rows();
                if(dim > 4 && counter[dim] < 10){

                    // Create filename in the output directory
                    std::string filename = std::to_string(dim) + "_" + std::to_string(counter[dim]) + ".scc";
                    std::filesystem::path specific_output_path = output_dir / filename;
                    
                    std::ofstream output_file(specific_output_path);
                    if (!output_file.is_open()) {
                        std::cerr << "Error: Unable to open output file " << specific_output_path << std::endl;
                        return;
                    } else {
                        ind.to_stream(output_file);
                        checked_list.emplace_back(current_index);
                        output_file.close();
                        counter[dim]++;
                        std::cout << "Large indecomposable Submodule at degree (" << degree.first << ", " << degree.second << ") computed and saved to: " << specific_output_path << std::endl;
                    }
                } 
            }
        }
    }
}

void large_induced_indecomp_submodules_from_sum(std::filesystem::path input_path, 
        std::filesystem::path output_dir) {
    std::vector<R2GradedSparseMatrix<int>> matrices;
    std::ifstream input_file(input_path);
    if (!input_file.is_open()) {
        std::cerr << "Error opening file: " << input_path << std::endl;
        return;
    }
    read_sccsum<int, std::ifstream>(matrices, input_file);
    vec<int> counter = vec<int>(50,0);
    for(auto& pres : matrices){
        large_induced_indecomp_submodules(pres, output_dir, counter);
    }
}

bool is_decomp_file(const std::filesystem::path& filepath) {
    return filepath.extension() == ".sccsum";
}

int main(int argc, char** argv) {
    std::filesystem::path input_path;

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <file_path> " << std::endl;
        input_path = "/home/wsljan/MP-Workspace/data/hypoxic_regions/hypoxic2_FoxP3_dim1_200x200_snapped.sccsum";
    } else {
        input_path = argv[1];
    }
    
    // Create output directory at the same location as input file
    std::filesystem::path output_dir = input_path.parent_path() / (input_path.stem().string() + "_induced");
    
    // Create the directory if it doesn't exist
    if (!std::filesystem::exists(output_dir)) {
        std::filesystem::create_directories(output_dir);
        std::cout << "Created output directory: " << output_dir << std::endl;
    }

    if (is_decomp_file(input_path)) {
        large_induced_indecomp_submodules_from_sum(input_path, output_dir);
    } else {
        R2GradedSparseMatrix<int> pres = R2GradedSparseMatrix<int>(input_path.string());
        large_induced_indecomp_submodules(pres, output_dir);
    }
    
    return 0;
}