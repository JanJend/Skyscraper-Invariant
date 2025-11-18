#include "filt_landscape.hpp"
#include <iostream>
#include <iomanip>

int main(int argc, char* argv[]) {
    if (argc > 6) {
        std::cerr << "Usage: " << argv[0] << " [<input.sky>] [<theta>] [<k>] \n";
        return 1;
    } 
    
    std::string input_file = (argc >= 2) ? argv[1] : "/home/wsljan/MP-Workspace/data/hypoxic_regions/hypoxic2_FoxP3_dim1_100x100_res_cut.sky";
    
    double theta = (argc >= 3) ? std::stod(argv[2]) : 0.0;
    int k = (argc >= 4) ? std::stoi(argv[3]) : 1;

    std::cout << "Computing landscape from file: " << input_file << " with theta = " << theta << " and k = " << k << std::endl;

    std::string output_file;
    if (argc >= 2) {
        // Remove extension and add "_landscape_theta_k.txt"
        std::ostringstream theta_stream;
        theta_stream << std::fixed << std::setprecision(2) << theta;
        std::string theta_str = theta_stream.str();
        
        size_t last_dot = input_file.find_last_of('.');
        output_file = (last_dot != std::string::npos) ? input_file.substr(0, last_dot) : input_file;
        output_file += "_landscape_" + theta_str + "_" + std::to_string(k) + ".txt";
    }
    
    
    
    try {
        hnf::GridData data = hnf::bars_from_sky(input_file);
        compute_landscape(data, output_file, theta, k);
        std::cout << "Landscape computed successfully\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}