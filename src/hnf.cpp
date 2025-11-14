#include "hnf.hpp"

namespace hnf {
// Uni_B1 template member function implementations
template<typename index>
Uni_B1<index>::Uni_B1(R2GradedSparseMatrix<index>&& d1_, bool is_minimal)
    : d1(std::move(d1_)) {
    if(is_minimal && d1.get_num_rows() == 1){
        assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
        d2 = R2GradedSparseMatrix<index>(d1.get_num_cols()-1, d1.get_num_cols());
        d2.data = vec< vec<index> >(d1.get_num_cols()-1);
        d2.row_degrees = d1.col_degrees;
        d2.col_degrees = vec<r2degree>(d1.get_num_cols()-1);
        r2degree last_degree = d1.col_degrees[0];
        for(index i = 1; i < d1.get_num_cols(); i++){
            r2degree join = Degree_traits<r2degree>::join(last_degree, d1.col_degrees[i]);
            d2.data[i] = {i -1, i};
            d2.col_degrees[i] = join;
            last_degree = d1.col_degrees[i];
        }
    } else {
        auto d1_copy = d1;
        d2 = d1_copy.graded_kernel();
    }
}

template<typename index>
Uni_B1<index>::Uni_B1(const R2GradedSparseMatrix<index>& d1_, bool is_minimal)
    : d1(d1_) {
    if(is_minimal && d1.get_num_rows() == 1){
        assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
        d2 = R2GradedSparseMatrix<index>(d1.get_num_cols()-1, d1.get_num_cols());
        d2.data = vec< vec<index> >(d1.get_num_cols()-1);
        d2.row_degrees = d1.col_degrees;
        d2.col_degrees = vec<r2degree>(d1.get_num_cols()-1);
        r2degree last_degree = d1.col_degrees[0];
        for(index i = 1; i < d1.get_num_cols(); i++){
            r2degree join = Degree_traits<r2degree>::join(last_degree, d1.col_degrees[i]);
            d2.data[i] = {i -1, i};
            d2.col_degrees[i] = join;
            last_degree = d1.col_degrees[i];
        }
    } else {
        auto d1_copy = d1;
        d2 = d1_copy.graded_kernel();
    }
}

template<typename index>
index Uni_B1<index>::dim_at(r2degree alpha) {
    index num_chains_0 = d1.num_rows_before(alpha);
    index num_chains_1 = d1.num_cols_before(alpha);
    index num_chains_2 = d2.num_cols_before(alpha);
    return num_chains_0 - num_chains_1 + num_chains_2;
}

template<typename index>
double Uni_B1<index>::area() const {
    auto [min, max] = d1.bounding_box();
    double base_area = d1.get_num_rows()*(max.first - min.first) * (max.second - min.second);
    for(const auto& degree : d1.col_degrees){
        base_area -= (max.first - degree.first) * (max.second - degree.second);
    }
    for(const auto& degree : d2.col_degrees){
        base_area += (max.first - degree.first) * (max.second - degree.second);
    }
    return base_area;
}

template<typename index>
double Uni_B1<index>::area(const r2degree& bound) const {
    auto [min, max] = d1.bounding_box();
    max = bound;
    double base_area = d1.get_num_rows()*(max.first - min.first) * (max.second - min.second);
    for(const auto& degree : d1.col_degrees){
        base_area -= (max.first - degree.first) * (max.second - degree.second);
    }
    for(const auto& degree : d2.col_degrees){
        base_area += (max.first - degree.first) * (max.second - degree.second);
    }
    return base_area;
}

template<typename index>
double Uni_B1<index>::area(const pair<r2degree>& bounds) const {
    auto [min, max] = d1.bounding_box();
    max = bounds.second;
    assert(max.first >= min.first);
    assert(max.second >= min.second);
    assert(min == d1.row_degrees[0]);
    r2degree range = bounds.second-bounds.first;
    double range_area = range.first * range.second;
    double base_area = d1.get_num_rows()*(max.first - min.first)*(max.second - min.second);
    for(const auto& degree : d1.col_degrees){
        assert(degree.first <= max.first);
        assert(degree.second <= max.second);
        assert(degree.first >= min.first);
        assert(degree.second >= min.second);
        base_area -= (max.first - degree.first) * (max.second - degree.second);
    }
    for(const auto& degree : d2.col_degrees){
        assert(degree.first <= max.first);
        assert(degree.second <= max.second);
        assert(degree.first >= min.first);
        assert(degree.second >= min.second);
        base_area += (max.first - degree.first) * (max.second - degree.second);
    }
    double normalised_area = base_area/range_area;
    return normalised_area;
}

template<typename index>
void Uni_B1<index>::compute_area_polynomial(const pair<r2degree>& bounds) {
    int num_rows = d1.get_num_rows();
    assert(std::is_sorted(d1.col_degrees.begin(), d1.col_degrees.end(), Degree_traits<r2degree>::lex_lambda()));
    bool test = true;

    r2degree range = bounds.second-bounds.first;
    double range_area = range.first * range.second;

    const r2degree& gen_degree = this->d1.row_degrees[0];
    area_polynomial.fill(0.0);

    r2degree gen_vector = bounds.second - gen_degree;
    area_polynomial[0] += num_rows * gen_vector.first * gen_vector.second;
    area_polynomial[1] -= num_rows * gen_vector.second;
    area_polynomial[2] -= num_rows * gen_vector.first;
    area_polynomial[3] += num_rows;

    for(const auto& degree : d1.col_degrees){
        assert( Degree_traits<r2degree>::smaller_equal(degree, bounds.second));
        assert( Degree_traits<r2degree>::smaller_equal(bounds.first, degree) );
        r2degree rel_vector = bounds.second - degree;
        area_polynomial[0] -= rel_vector.first* rel_vector.second;
    }

    for(const auto& degree :  d1.col_degrees){
        if(degree.first == gen_degree.first){
            assert(bounds.second.second - degree.second >= 0);
            area_polynomial[1] += (bounds.second.second - degree.second);
            if( area_polynomial[1] > 0){
                std::cout << "Area polynomial[1] is positive: " << area_polynomial[1] << std::endl;
                this->d1.print_graded();
                std::cout << "Bounds: " << bounds.first.first << " " << bounds.first.second << " to "
                    << bounds.second.first << " " << bounds.second.second << std::endl;
                std::cout << "Range area: " << range_area << std::endl;
                assert(false);
            }
        }
        if(degree.second == gen_degree.second){
            assert(bounds.second.first - degree.first >= 0);
            area_polynomial[2] += (bounds.second.first - degree.first);
            if( test && area_polynomial[2] > 0){
                std::cout << "Area polynomial[2] is positive: " << area_polynomial[1] << std::endl;
                this->d1.print_graded();
                std::cout << "Bounds: " << bounds.first.first << " " << bounds.first.second << " to "
                    << bounds.second.first << " " << bounds.second.second << std::endl;
                std::cout << "Range area: " << range_area << std::endl;
                assert(false);
            }
        }
    }

    for(const auto& degree : d2.col_degrees){
        assert( Degree_traits<r2degree>::smaller_equal(degree, bounds.second));
        assert( Degree_traits<r2degree>::smaller_equal(bounds.first, degree) );
        r2degree rel_vector = bounds.second - degree;
        area_polynomial[0] += rel_vector.first* rel_vector.second;
    }
    for(index i = 0; i < 4; i++){
        area_polynomial[i] /= range_area;
    }
}

template<typename index>
double Uni_B1<index>::evaluate_area_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    return area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + area_polynomial[3]*x*y;
}

template<typename index>
double Uni_B1<index>::evaluate_slope_polynomial(r2degree d) {
    double& x = d.first;
    double& y = d.second;
    return 1/ (area_polynomial[0] + area_polynomial[1]*x + area_polynomial[2]*y + area_polynomial[3]*x*y);
}

template<typename index>
double Uni_B1<index>::slope() const {
    double area = this->area();
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. Consider passing a bound." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}

template<typename index>
double Uni_B1<index>::slope(const r2degree& bound) const {
    double area = this->area(bound);
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. The bound you passed is insufficient." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}

template<typename index>
double Uni_B1<index>::slope(const pair<r2degree>& bounds) const {
    double area = this->area(bounds);
    if(area == 0){
        std::cerr << "Area is zero, slope will be infinite. The bound you passed is insufficient." << std::endl;
    }
    double slope_value = (static_cast<double>(d1.get_num_rows())/area);
    return slope_value;
}


// slope_comparator
bool slope_comparator::operator()(const Module_w_slope& X, const Module_w_slope& Y) const noexcept {
    return X.second > Y.second;
}

// Free functions (non-template)
Module_w_slope find_scss_bruteforce(const R2Mat& X,
        vec<vec<SparseMatrix<int>>>& subspaces,
        R2Mat& max_subspace,
        const pair<r2degree>& bounds) {
    int k = X.get_num_rows();
    Uni_B1<int> scss = Uni_B1<int>(X);
    double max_slope = scss.slope(bounds);
    if(k == 1){
        // Nothing to do?
    } else {
        assert(k < 6);
        if(subspaces.size() < k){
            std::cerr << "Have not loaded enough subspaces" << std::endl;
            std::exit(1);
        }
        for(auto ungraded_subspace : subspaces[k-1]){
            int num_gens = ungraded_subspace.get_num_cols();
            R2Mat subspace = R2Mat(ungraded_subspace);
            subspace.row_degrees = X.row_degrees;
            subspace.col_degrees = vec<r2degree>(num_gens, X.row_degrees[0]);
            assert(subspace.get_num_rows() == X.get_num_rows());
            assert(subspace.get_num_cols() == num_gens);
            R2Mat submodule_pres = X.submodule_generated_by(subspace);
            Uni_B1<int> res(submodule_pres);
            double slope = res.slope(bounds);
            if(slope > max_slope){
                max_slope = slope;
                scss = std::move(res);
                max_subspace = subspace;
            }
        }
    }
    return std::make_pair(scss, max_slope);
}

HN_factors skyscraper_invariant(Block_list& summands,
        vec<vec<SparseMatrix<int>>>& subspaces,
        const pair<r2degree>& bounds) {
    HN_factors result;
    for(Block X : summands){
        if(X.get_num_rows() >= 4){
            std::cout << "  Warning: Computing HNF for a module of dimension "
                << X.get_num_rows() << std::endl;
        }

        while(X.get_num_rows() > 1){
            R2Mat max_subspace;
            result.emplace_back(find_scss_bruteforce(X, subspaces, max_subspace, bounds));
            if(result.back().first.d1.get_num_rows() == X.get_num_rows()){
                break;
            } else {
                X.quotient_by(max_subspace);
            }
        }
        if(X.get_num_rows() == 1){
            Uni_B1<int> res(X);
            double slope = res.slope(bounds);
            result.emplace_back(std::make_pair(res, slope));
        } else if (X.get_num_rows() == 0){
            // Nothing to do.
        }
    }
    return result;
}

void skyscraper_invariant(const R2Mat& input,
    HN_factors& result,
    vec<vec<SparseMatrix<int>>>& subspaces,
    const pair<r2degree>& bounds) {
    
    R2Mat X = input;

    if(X.get_num_rows() > 1){
        
       // X.to_stream_r2(std::cout);
    }

    if(X.get_num_rows() >= 6){
        std::cout << "  Warning: Computing HNF for a module of dimension "
            << X.get_num_rows() << std::endl;
    }
    while(X.get_num_rows() > 1){
        R2Mat subspace;
        result.emplace_back(find_scss_bruteforce(X, subspaces, subspace, bounds));
        if(result.back().first.d1.get_num_rows() == X.get_num_rows()){
            break;
        } else {
            X.quotient_by(subspace);
        }
    }
    if(X.get_num_rows() == 1){
        Uni_B1<int> res(X);
        double slope = res.slope(bounds);
        result.emplace_back(std::make_pair(res, slope));
    } else if (X.get_num_rows() == 0){
        // Nothing to do.
    }
}

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
        slopes.push_back(hn_factor.second);
        test_slopes.push_back(it->second);
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
            composition_factors[l].first.d1.print_graded();
            std::cout << "  Test_summand: " << std::endl;
            test_factors[l].first.d1.print_graded();
            assert(false);
        }
    }
}

// Dynamic_HNF
Dynamic_HNF::Dynamic_HNF() {
    indecomposable_summands = vec<vec<Uni_B1<int>>>();
    grid_ind_dimensions = vec<int>();
}

void Dynamic_HNF::compute_HNF_row(aida::AIDA_functor& decomposer,
        R2Mat& M, int& y_index, pair<r2degree> slope_bounds) {
    assert(y_index > -1);
    double y_coordinate = M.y_grid[y_index];
    int x_length = M.x_grid.size();

    indecomposable_summands = vec<vec<Uni_B1<int>>>(x_length, vec<Uni_B1<int>>());
    int max_dim = 0;

    for(int x_index = 0; x_index < x_length; x_index++){
        r2degree grid_point = {M.x_grid[x_index], y_coordinate};
        auto M_induced = M.submodule_generated_at(grid_point);
        if(M_induced.get_num_rows() == 0){
            // Nothing to do?
        } else if (M_induced.get_num_rows() == 1){
            Uni_B1<int> res(std::move(M_induced));
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
                indecomposable_summands[x_index].emplace_back(Uni_B1<int>(std::move(sub_M)));
                indecomposable_summands[x_index].back().compute_area_polynomial(slope_bounds);
                grid_ind_dimensions.push_back(sub_M.get_num_rows());
            }
        }
    }
    if (max_dim >= 6){
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
