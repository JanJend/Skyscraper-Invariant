#include "hnf.hpp"

namespace hnf {


void calculate_stats(const std::vector<int>& all_dimensions) {
    if (all_dimensions.empty()) {
        std::cout << "The vector is empty!" << std::endl;
        return;
    }

    int max_value = *std::max_element(all_dimensions.begin(), all_dimensions.end());

    double sum = std::accumulate(all_dimensions.begin(), all_dimensions.end(), 0);
    double average = sum / all_dimensions.size();

    double squared_diff_sum = 0;
    for (int val : all_dimensions) {
        squared_diff_sum += (val - average) * (val - average);
    }
    double variance = squared_diff_sum / all_dimensions.size();
    double standard_deviation = std::sqrt(variance);
    int count_not_one = std::count_if(all_dimensions.begin(), all_dimensions.end(), 
                                   [](int val) { return val != 1; });
    double percentage_not_one = (static_cast<double>(count_not_one) / all_dimensions.size()) * 100.0;

    std::cout << std::fixed << std::setprecision(8);
    std::cout << "Percentage not 1: " << percentage_not_one << "%" << std::endl;
    std::cout << "Maximum: " << max_value << std::endl;
    std::cout << "Average: " << average << std::endl;
    std::cout << "Standard Deviation: " << standard_deviation << std::endl;
}

vec<r2degree> get_grid_diagonal(pair<r2degree> bounds, int grid_length) {
    vec<r2degree> grid_diagonal = vec<r2degree>();
    double x_min = bounds.first.first;
    double x_max = bounds.second.first;
    double y_min = bounds.first.second;
    double y_max = bounds.second.second;

    double x_step = (x_max - x_min) / (grid_length - 1);
    double y_step = (y_max - y_min) / (grid_length - 1);

    for (int i = 0; i < grid_length; ++i) {
        grid_diagonal.push_back({x_min + i * x_step, y_min + i * y_step});
    }
    return grid_diagonal;
}

r2degree get_grid_step(const r2degree& lower_bound, const r2degree& upper_bound,
    const int& grid_length_x, const int& grid_length_y) {
    const double& x_min = lower_bound.first;
    const double& x_max = upper_bound.first;
    const double& y_min = lower_bound.second;
    const double& y_max = upper_bound.second;

    double x_step = (x_max - x_min) / (grid_length_x - 1);
    double y_step = (y_max - y_min) / (grid_length_y - 1);

    return {x_step, y_step};
}

void write_slopes_to_csv(const vec<vec<double>>& slopes,
        const vec<r2degree>& grid_points,
        const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file");
    }

    for (size_t i = 0; i < slopes.size(); ++i) {
        file << grid_points[i] << ";";

        const auto& slopes_at_degree = slopes[i];
        for (size_t j = 0; j < slopes_at_degree.size(); ++j) {
            file << slopes_at_degree[j];
            if (j < slopes_at_degree.size() - 1) {
                file << ",";
            }
        }    
        file << "\n";
    }
    file.close();
}

void show_progress_bar(int& i, int& total, std::string& name) {
    static int last_percent = -1;
    int percent = static_cast<int>(static_cast<double>(i+1) / total * 100);
    if (percent != last_percent) {
        int num_symbols = percent / 2;
        std::cout << "\r" << i + 1 << " " << name << " : [";
        for (int j = 0; j < 50; ++j) {
            if (j < num_symbols) {
                std::cout << "#";
            } else {
                std::cout << " ";
            }
        }
        std::cout << "] " << percent << "%";
        std::flush(std::cout);
        last_percent = percent;
    }
    if (i >= total) {
        std::cout << std::endl;
    }
}

bool essentially_equal(double a, double b, double relTol, double absTol) {
    return std::fabs(a - b) <= std::max(relTol * std::max(std::fabs(a), std::fabs(b)), absTol);
}

void compare_slopes_test(r2degree current_grid_degree,
    r2degree local_grid_degree,
    const HN_factors& composition_factors,
    const HN_factors& test_factors,
    int i, int j, int k) {
    vec<double> slopes;
    vec<double> test_slopes;
    assert(test_factors.size() == composition_factors.size());
    auto it = test_factors.begin();
    for(auto& hn_factor : composition_factors){
        slopes.push_back(hn_factor.slope_value);
        test_slopes.push_back(it->slope_value);
        it++;
    }
    for(size_t l = 0; l < slopes.size(); l++){
        if(!essentially_equal(slopes[l], test_slopes[l], 1e-6, 1e-8)){
            std::cout << std::fixed << std::setprecision(12);
            std::cout << "  Slope mismatch in compare_slopes: " << slopes[l] << " vs. " << test_slopes[l] << std::endl;
            std::cout << "  Difference: " << slopes[l] - test_slopes[l] << std::endl;
            std::cout << "  Current grid degree: " << current_grid_degree << std::endl;
            std::cout << "  Local grid degree: " << local_grid_degree << std::endl;
            std::cout << "  i: " << i << ", j: " << j << ", k: " << k << std::endl;
            std::cout << "  Summand: " << std::endl;
            composition_factors[l].d1.print_graded();
            std::cout << "  Test_summand: " << std::endl;
            test_factors[l].d1.print_graded();
            assert(false);
        }
    }
}

// Dynamic_HNF
Dynamic_HNF::Dynamic_HNF() {
    indecomposable_summands = vec<vec<Uni_B1>>();
    grid_ind_dimensions = vec<int>();
}

void Dynamic_HNF::compute_HNF_row(aida::AIDA_functor& decomposer,
        R2Mat& M, int& y_index, pair<r2degree> slope_bounds,
        const vec<vec<vec<SparseMatrix<int>>>>& subspaces) {
    assert(y_index > -1);
    double y_coordinate = M.y_grid[y_index];
    int x_length = M.x_grid.size();
    double y_next;
    if(y_index < M.y_grid.size()-1){
        y_next = M.y_grid[y_index+1];
    } else {
        y_next = slope_bounds.second.second;
    }
    indecomposable_summands = vec<vec<Uni_B1>>(x_length, vec<Uni_B1>());
    int max_dim = 0;

    for(int x_index = 0; x_index < x_length; x_index++){
        r2degree grid_point = {M.x_grid[x_index], y_coordinate};
        double next_x;
        if(x_index < M.x_grid.size()-1){
            next_x = M.x_grid[x_index+1];
        } else {
            next_x = slope_bounds.second.first;
        }
        auto M_induced = M.submodule_generated_at(grid_point);
        if(M_induced.get_num_rows() == 0){
            // Nothing to do?
        } else if (M_induced.get_num_rows() == 1){
            Uni_B1 res(std::move(M_induced));
            indecomposable_summands[x_index].push_back(res);
            indecomposable_summands[x_index].back().compute_area_polynomial(slope_bounds);
            grid_ind_dimensions.push_back(1);
        } else {
            aida::Block_list sub_M_list;
            M_induced.compute_col_batches();
            decomposer(M_induced, sub_M_list);
            for(Block sub_M : sub_M_list){
                if(sub_M.get_num_rows() > max_dim){
                    max_dim = sub_M.get_num_rows();
                }
                indecomposable_summands[x_index].emplace_back(Uni_B1(std::move(sub_M)));
                Uni_B1& current_summand =  indecomposable_summands[x_index].back();
                current_summand.compute_area_polynomial(slope_bounds);
                grid_ind_dimensions.push_back(sub_M.get_num_rows());
                int dim = current_summand.d1.get_num_rows();
                if(false){
                    if(dim > 2){
                        current_summand.d1.to_stream_r2(std::cout);
                    }
                }

                if(false){
                    r2degree upper_grid_corner = {next_x, y_next};
                    // current_summand.compute_slope_subdivision(slope_bounds, subspaces, grid_point, upper_grid_corner);
                }
            }
        }
    }
    if (max_dim >= 7){
        // std::cout << " Careful, there are high-dimensional summands which might slow down HNF computation excessively." << std::endl;
    }  
}






// All of the following should normally go in a main.cpp file.

namespace fs = std::filesystem;

void display_help() {
    std::cout << "Usage: ./aida <input_file> [options]\n"
              << "Options:\n"
              << "  -h, --help           Display this help message\n"
              << "  -g, --diagonal       Save a copy where each subquotient is restricted to the diagonal to compute landscapes\n"
              << "  -d, --is_decomposed  Specify if the input is already decomposed\n"
              << "  -v, --version        Display version information\n"
              << "  -b, --bruteforce     Stops hom-space calculation and thus most optimisation. \n"
              << "  -s, --sort           Lexicographically sorts the relations of the input\n"
              << "  -e, --exhaustive     Always iterates over all decompositions of a batch\n"
              << "  -t, --statistics     Show statistics about indecomposable summands\n"
              << "  -r, --runtime        Show runtime statistics and timers\n"
              << "  -p, --progress       Turn off progressbar\n"
              << "  -c, --basechange     Save base change\n"
              << "  -o, --output <file>  Specify output file\n"
              << "  -l, --less_console   Suppreses most console output\n"
              << "  -m, --compare_b      Compares with -b at runtime, then runs with only -b and compares.\n"
              << "  -a, --compare_e      Compares exhaustive and brute force at runtime.\n"
              << "  -i, --compare_hom    Compares optimised and non-opt hom space calculation at runtime.\n"
              << "  -j, --no_hom_opt     Does not use the optimised hom space calculation.\n"
              << "  -w, --no_col_sweep   Does not use the column sweep optimisation.\n"
              << "  -f, --alpha       Turns the computation of alpha-homs on.\n"
              << "  -x, -test_files          Runs the algorithm on some test files.\n"
              << "      <file> is optional and will default to the <input_file> with _decomposed appended\n"
              << "      You can pass relative and absolute paths as well as only a directory."
              << "Further Instructions: \n Make sure that the inputfile is a (sequence of) scc or firep presentations that are minimised.\n"
              << std::endl;
}

void display_version() {
    std::cout << "Skyscraper Invariant version 0.1 -- 7th Oct 2025\n";
}



void write_to_file(std::ostringstream& ostream, std::string& output_file_path, std::string& input_directory, std::string& file_without_extension, std::string& extension, std::string& output_string){

    if(output_string.empty()){
        output_file_path = input_directory + "/" + file_without_extension + extension;
    } else {
        std::filesystem::path output_path(output_string);
        if (output_path == ".") {
            output_file_path = std::filesystem::current_path().string() + "/" + file_without_extension +  extension;
        } else if (output_path.is_relative()) {
            output_file_path = std::filesystem::current_path().string() + "/" + output_string;
        } else if (std::filesystem::is_directory(output_path)) {
            output_file_path = output_path.string() + "/" + file_without_extension + extension;
        } else if (output_path.is_absolute()) {
            output_file_path = output_string;
        } else {
            output_file_path = input_directory + "/" + output_string;
        }
    }

    std::filesystem::create_directories(std::filesystem::path(output_file_path).parent_path());

    std::ofstream file_out(output_file_path);
    if(file_out.is_open()){
        file_out << ostream.str();
        file_out.close();
        std::cout << "HN filtration written to " << output_file_path << std::endl;
    } else {
        std::cout << "Error: Could not write HN filtration to file: " << output_file_path << std::endl;
    }
}

} // namespace hnf
