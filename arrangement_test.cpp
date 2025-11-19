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
    Uni_B1<int> M_res(M);
    auto bounding_box = M_res.d1.bounding_box();
    r2degree cell_boundary = bounding_box.second;
    // Notice that right now, this is the whole support, and not just the grid cell.
    auto supspaces = all_sparse_subspaces(k);
    M_res.compute_slope_subdivision(bounding_box, supspaces, cell_boundary);
    M_res.slope_subdiv->export_to_svg(output_path.string());
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