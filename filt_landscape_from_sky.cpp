#include "filt_landscape.hpp"
#include <iostream>
#include <iomanip>

using namespace hnf;

int main(int argc, char* argv[]) {
    if (argc > 6 || argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.sky> <theta> <k> <diff> <theta_prime> \n";
        std::cerr << "  <input.sky> : Path to the input skyscraper file.\n";
        std::cerr << "  <theta>     : double (Optional) Filtration parameter (default: 0.0).\n";
        std::cerr << "  <k>         : int (Optional) Landscape level (default: 1).\n";
        std::cerr << "  <diff>      : bool (Optional) 'true' to compute difference landscape.\n";
        std::cerr << "  <theta_prime> : double (Optional) Second filtration parameter for difference landscape (default: 0.0).\n";
        
        return 1;
    } 
    
    std::string input_file = (argc >= 2) ? argv[1] : "/home/wsljan/MP-Workspace/data/hypoxic_regions/hypoxic2_FoxP3_dim1_200x200_res.sky";
    
    double theta = (argc >= 3) ? std::stod(argv[2]) : 0.0;
    int k = (argc >= 4) ? std::stoi(argv[3]) : 1;
    bool diff = (argc >= 5) ? (std::string(argv[4]) == "true") : false;
    double theta_prime = (argc >= 6) ? std::stod(argv[5]) : 0.0;
    std::cout << "Computing" << (diff ? " difference" : "") << " landscape from file: " << input_file << " with theta = " << theta << ", theta_prime = " << theta_prime << " and k = " << k << std::endl;

    std::string output_file;
    if (argc >= 2) {
        // Remove extension and add "_landscape_theta_k.txt"
        std::ostringstream theta_stream;
        theta_stream << std::fixed << std::setprecision(2) << theta;
        std::string theta_str = theta_stream.str();
        std::ostringstream theta_prime_stream;
        theta_prime_stream << std::fixed << std::setprecision(2) << theta_prime;
        std::string theta_prime_str = theta_prime_stream.str();
        size_t last_dot = input_file.find_last_of('.');
        output_file = (last_dot != std::string::npos) ? input_file.substr(0, last_dot) : input_file;
        output_file += "_landscape_" + theta_str + (diff ? "_diff" + theta_prime_str : "") + ".png";
    }
    
    try {
        hnf::GridData data = hnf::bars_from_sky(input_file);
        std::vector<std::vector<std::vector<double>>> landscape;
        if(diff){
            landscape = hnf::compute_difference_landscape(data, theta, theta_prime, k);
            std::cout << "Difference landscape computed successfully\n";
        } else {
            landscape = hnf::compute_landscape(data, theta, k);
            std::cout << "Landscape computed successfully\n";
        }
        hnf::write_landscape_png(landscape, output_file);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}