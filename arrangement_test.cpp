#include "subdivision.hpp"

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

    R2PModule module(input_path.string());
    R2Mat& M = module.mutable_presentation();
    int k = M.get_num_rows();
    pair<r2degree> bounding_box = geometry_bounds(M);

    // enlarging bounding box to see whole arrangement
    
    assert(M.grid(0).size() >= 2 && M.grid(1).size() >= 2);
    double cell_end_x = M.grid(0)[1];
    double cell_end_y = M.grid(1)[1];
    r2degree cell_start = {M.grid(0)[0], M.grid(1)[0]};
    r2degree cell_boundary = {cell_end_x, cell_end_y};
    cell_boundary = bounding_box.second;
    r2degree range = bounding_box.second - bounding_box.first;
    r2degree padding = {0.2,0.2};
    r2degree large_padding = {1.2, 1.2};
    cell_start = cell_start - large_padding;
    cell_boundary = cell_boundary + large_padding;
    // cell_start = {-2, -2};
    // cell_boundary = {2, 2};
    std::cout << "Computing arrangement in the box: (" << cell_start[0] << ", " << cell_start[1]
        << ") to (" << cell_boundary[0] << ", " << cell_boundary[1] << ")\n";
    Uni_B1 M_res(module);
    auto subspaces = all_sparse_subspaces(k);
    Slope_subdivision slope_subdiv = compute_slope_subdivision(M_res, bounding_box, subspaces, cell_start, cell_boundary);
    slope_subdiv.export_to_svg(output_path.string(), bounding_box.first[0], bounding_box.first[1]);
}

int main(int argc, char** argv) {
    std::string filepath;
    int optional_value = 0;  // default value
    std::string extension = ".svg";
    std::filesystem::path output_path;
    std::string suffix;
    std::filesystem::path input_path;

    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <file_path> [optional_value] \n";
        std::cerr << "  optional_value = 0 \n";
        return 1;
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
