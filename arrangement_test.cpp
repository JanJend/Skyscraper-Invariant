#include "hnf_uni_b1.hpp"

namespace fs = std::filesystem;

using namespace graded_linalg;
using namespace hnf;

void get_arrangement(std::filesystem::path& input_path,
    std::filesystem::path& output_path,
    int optional_value) {
    
    std::ifstream istream(input_path);
    if (!istream.is_open()) {
        std::cerr << "Error: Could not open input file " << input_path << std::endl;
        return;
    }

    R2GradedSparseMatrix<int> M(input_path.string());
    int k = M.get_num_rows();
    auto perm = M.compute_grid_representation();
    pair<r2degree> bounding_box = M.bounding_box();

    // enlarging bounding box to see whole arrangement
    

    assert(M.x_grid.size() >= 2 && M.y_grid.size() >= 2);
    double cell_end_x = M.x_grid[1];
    double cell_end_y = M.y_grid[1];
    r2degree cell_start = {M.x_grid[0], M.y_grid[0]};
    r2degree cell_boundary = {cell_end_x, cell_end_y};
    cell_boundary = bounding_box.second;
    r2degree range = bounding_box.second - bounding_box.first;
    cell_start = cell_start - 0.2*range;
    // cell_start = {-2, -2};
    // cell_boundary = {2, 2};
    std::cout << "Computing arrangement in the box: (" << cell_start.first << ", " << cell_start.second 
        << ") to (" << cell_boundary.first << ", " << cell_boundary.second << ")\n";
    Uni_B1 M_res(M);
    auto subspaces = all_sparse_subspaces(k);
    M_res.compute_slope_subdivision(bounding_box, subspaces, cell_start, cell_boundary);
    M_res.slope_subdiv->export_to_svg(output_path.string(), bounding_box.first.first, bounding_box.first.second);
}

int main(int argc, char** argv) {
    std::string filepath;
    int optional_value = 0;  // default value
    std::string extension = ".svg";
    std::filesystem::path output_path;
    std::string suffix;
    std::filesystem::path input_path;

    if (argc < 2 || argc > 5) {
        std::cerr << "Usage: " << argv[0] << " <file_path> [optional_value] \n";
        std::cerr << "  optional_value = 0 \n";
        filepath = "/home/wsljan/MP-Workspace/Skyscraper-Invariant/example_files/indecomps_at/stable_example_paper.scc";
        input_path = std::filesystem::path(filepath);
    } else {
        filepath = argv[1];
        input_path = std::filesystem::path(filepath);
    }
    if (argc >= 3) {
        try {
            optional_value = std::stoi(argv[2]);
        } catch (const std::exception& e) {
            std::cerr << "Error: Invalid integer argument" << std::endl;
            return 1;
        }
    }
    
    suffix = "_" + std::to_string(optional_value);
    std::string modified_path = insert_suffix_before_extension(filepath, suffix, extension);
    output_path = std::filesystem::path(modified_path);
    
    get_arrangement(input_path, output_path, optional_value);

    return 0;
}