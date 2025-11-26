#include <grlina/graded_linalg.hpp>
#include <iostream>
#include <filesystem>

using namespace graded_linalg;


void compute_random_uni_b1(std::filesystem::path output_path, const int optional_value = 10){
    using R2Matrix = R2GradedSparseMatrix<int>;
    int dim = optional_value;
    int num_rel = std::pow(2, dim+3);
    int extra_fill_chance = std::pow(2, dim);
    int perc = 100/dim;
    std::cout << "Generating random Uni B1 presentation with dimension " << dim 
        << ", " << num_rel << " relations and fill percentage " << perc << std::endl;
    R2Matrix M(num_rel, dim, "Random", perc);
    M.row_degrees = vec<r2degree>(M.get_num_rows(), {0.0, 0.0});
    // Set random col degrees in the box [0,1] x [0,1]
    for(int i = 0; i < M.get_num_cols(); i++){
        M.col_degrees[i] = {static_cast<double>(rand()%100)/100.0, static_cast<double>(rand()%100)/100.0};
        double fill_perc = static_cast<double>(rand()%100)/100.0;
        if( fill_perc < 1/static_cast<double>(extra_fill_chance) ){
            double fill_amount = static_cast<double>(rand()%100)/100.0;
            for(int j = 0; j < dim; j++){
                if(static_cast<double>(rand()%100)/100.0 < fill_amount){
                    M.data[i].push_back(j);
                }
            }
            std::sort(M.data[i].begin(), M.data[i].end());
            convert_mod_2(M.data[i]);
        }
    }
    // Make the module bounded
    for(int i = 0; i < dim ; i++){
        M.data.push_back( {i} );
        M.data.push_back( {i} );
        M.col_degrees.push_back( {1.0, 0.0} );
        M.col_degrees.push_back( {0.0, 1.0} );
    }
    M.set_num_cols(M.data.size());

    M.print_graded();
    M.sort_columns_lexicographically();
    M.minimize();
    M.print_graded();
    std::ofstream output_file(output_path);
    if (!output_file.is_open()) {
        std::cerr << "Error: Unable to open output file " << output_path << std::endl;
        return;
    } else {
        // M.to_stream_r2(output_file);
        output_file.close();
        std::cout << "Random Uni B1 presentation computed and saved to: " << output_path << std::endl;
    }
}

int main(int argc, char** argv) {

    int optional_value = 3;  // default value
    std::filesystem::path output_path;
    std::string suffix;

    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <dimension> [output_path] \n";
    } else {
        try {
            optional_value = std::stoi(argv[1]);
            if(optional_value <= 0){
                throw std::invalid_argument("dimension must be positive");
            }
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid integer argument" << std::endl;
            return 1;
        }
    }
    
    if(argc >= 3) {
        try {
            output_path = std::filesystem::path(argv[2]);
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid output path argument" << std::endl;
            return 1;
        }
    } else {
        std::string extension = ".scc";
        if (optional_value != 0){
            suffix = "dim_" + std::to_string(optional_value);
        } else {
            std::cerr << "Error: Empty output." << std::endl;
        }   
        std::string modified_path = "random_dim_" + std::to_string(optional_value) + extension;
        output_path = std::filesystem::path(modified_path);
    }
    
    compute_random_uni_b1(output_path, optional_value);

    return 0;
}