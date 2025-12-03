#include "hnf_interface.hpp"

namespace fs = std::filesystem;

struct ProgramConfig {
    aida::AIDA_functor decomposer;
    bool write_output = false;
    bool diagonal_output = false;
    bool show_indecomp_statistics = false;
    bool show_runtime_statistics = false;
    bool test_files = false;
    bool is_decomposed = false;
    bool dynamic_grid = true;
    bool subdivision = false;
    int grid_length_x = 200;
    int grid_length_y = 200;
    int grassmann_value = -1;
    std::string output_string;
};

struct FileInfo {
    std::string matrix_path;
    std::string input_directory;
    std::string filename;
    std::string file_without_extension;
    std::string extension = ".sky";
};

void initialize_decomposer_config(aida::AIDA_functor& decomposer) {
    decomposer.config.brute_force = false; // There is currently a bug in aida when this is false.
    decomposer.config.exhaustive = false;
    decomposer.config.sort = false;
    decomposer.config.sort_output = true;
    decomposer.config.alpha_hom = false;
    decomposer.config.progress = true;
    decomposer.config.show_info = true;
    decomposer.config.compare_both = false;
    decomposer.config.exhaustive_test = false;
    decomposer.config.save_base_change = false;
    decomposer.config.turn_off_hom_optimisation = false;
}

bool handle_interactive_input(int& argc, char**& argv) {
    std::cerr << "No input file specified. Please provide an input file." << std::endl;
    hnf::display_help();
    std::cout << "Please provide options/arguments: ";
    
    std::string input;
    std::getline(std::cin, input);
    
    std::vector<std::string> args;
    args.push_back(argv[0]);
    std::istringstream iss(input);
    std::string token;
    while (iss >> token) {
        args.push_back(token);
    }
    
    argc = args.size();
    argv = new char*[argc + 1];
    for (size_t i = 0; i < args.size(); ++i) {
        argv[i] = new char[args[i].size() + 1];
        std::strcpy(argv[i], args[i].c_str());
    }
    argv[argc] = nullptr;
    optind = 1;
    
    return true;
}

bool parse_resolution(const std::string& res_arg, int& grid_x, int& grid_y) {
    size_t comma_pos = res_arg.find(',');
    if (comma_pos == std::string::npos) {
        std::cerr << "Error: Resolution argument must be in the format 'x,y'." << std::endl;
        return false;
    }
    grid_x = std::stoi(res_arg.substr(0, comma_pos));
    grid_y = std::stoi(res_arg.substr(comma_pos + 1));
    return true;
}

bool parse_command_line(int argc, char** argv, ProgramConfig& config) {
    static struct option long_options[] = {
        {"help", no_argument, 0, 'h'},
        {"version", no_argument, 0, 'v'},
        {"output", optional_argument, 0, 'o'},
        {"exhaustive", no_argument, 0, 'e'},
        {"statistics", no_argument, 0, 's'},
        {"runtime", no_argument, 0, 't'},
        {"progress", no_argument, 0, 'p'},
        {"basechange", no_argument, 0, 'c'},
        {"less_console", no_argument, 0, 'l'},
        {"no_hom_opt", no_argument, 0, 'j'},
        {"alpha", no_argument, 0, 'f'},
        {"test_files", no_argument, 0, 'x'},
        {"is_decomposed", no_argument, 0, 'd'},
        {"diagonal", no_argument, 0, 'g'},
        {"resolution", required_argument, 0, 'r'},
        {"dynamic_grid", no_argument, 0, 'y'},
        {"subdivision", no_argument, 0, 'u'},
        {"grassmann", required_argument, 0, 'k'},
        {0, 0, 0, 0}
    };
    
    int opt;
    int option_index = 0;
    
    while ((opt = getopt_long(argc, argv, "ho::gestr:pclfjxdyk:u", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'h':
                hnf::display_help();
                return false;
            case 'v':
                hnf::display_version();
                return false;
            case 'o':
                config.write_output = true;
                if (optarg) {
                    config.output_string = std::string(optarg);
                } else if (optind < argc && argv[optind][0] != '-') {
                    config.output_string = std::string(argv[optind]);
                    optind++;
                } else {
                    config.output_string.clear();
                }
                break;
            case 's':
                config.show_indecomp_statistics = true;
                break;
            case 'e':
                config.decomposer.config.exhaustive = true;
                break;
            case 't':
                config.show_runtime_statistics = true;
                break;
            case 'r':
                if (!optarg) {
                    std::cerr << "Error: No resolution argument provided." << std::endl;
                    return false;
                }
                if (!parse_resolution(optarg, config.grid_length_x, config.grid_length_y)) {
                    return false;
                }
                break;
            case 'p':
                config.decomposer.config.progress = false;
                break;
            case 'c':
                config.decomposer.config.save_base_change = true;
                break;
            case 'l':
                config.decomposer.config.show_info = false;
                break;
            case 'j':
                config.decomposer.config.turn_off_hom_optimisation = true;
                break;
            case 'f':
                config.decomposer.config.alpha_hom = true;
                break;
            case 'x':
                config.test_files = true;
                break;
            case 'd':
                config.is_decomposed = true;
                std::cout << "Input file is assumed to be already decomposed." << std::endl;
                break;
            case 'g':
                config.diagonal_output = true;
                break;
            case 'y':
                config.dynamic_grid = false;
                break;
            case 'u':
                config.subdivision = true;
                break;
            case 'k':
                if (!optarg) {
                    std::cerr << "Error: --grassmann requires an integer argument." << std::endl;
                    return false;
                }
                config.grassmann_value = std::stoi(optarg);
                break;
            default:
                return false;
        }
    }
    
    return true;
}

FileInfo resolve_input_file(int argc, char** argv, bool test_files, bool& is_decomposed) {
    FileInfo file_info;
    
    if (optind < argc) {
        std::filesystem::path fs_path(argv[optind]);
        file_info.matrix_path = fs_path.is_relative() 
            ? std::filesystem::current_path().string() + "/" + argv[optind]
            : argv[optind];
        file_info.input_directory = fs_path.parent_path().string();
        file_info.filename = fs_path.filename().string();
        
        size_t dot_position = file_info.filename.find_last_of('.');
        file_info.file_without_extension = (dot_position == std::string::npos) 
            ? file_info.filename 
            : file_info.filename.substr(0, dot_position);
    } else if (test_files) {
        // Do nothing - test files handled separately
    } else {
        // No file provided and not in test mode - use default test file
        fs::path default_test_file = "/home/wsljan/MP-Workspace/data/hypoxic_regions/hypoxic2_FoxP3_dim1_200x200_res.sccsum";
        
        file_info.matrix_path = default_test_file.string();
        file_info.input_directory = default_test_file.parent_path().string();
        file_info.filename = default_test_file.filename().string();
        
        size_t dot_position = file_info.filename.find_last_of('.');
        file_info.file_without_extension = (dot_position == std::string::npos) 
            ? file_info.filename 
            : file_info.filename.substr(0, dot_position);
        
        is_decomposed = true;
        std::cout << "No input file specified. Running on test file: " << file_info.matrix_path << std::endl;
    }
    
    return file_info;
}

bool process_input_file(const FileInfo& file_info, ProgramConfig& config, std::ostringstream& ostream) {
    std::ifstream istream(file_info.matrix_path);
    if (!istream.is_open()) {
        std::cerr << "Error: Could not open input file: " << file_info.matrix_path << std::endl;
        return false;
    }
    
    std::cout << (config.is_decomposed 
        ? "Running HNF on already decomposed input file: " 
        : "First decomposing with AIDA.") + file_info.filename << std::endl;
    
    ostream << std::fixed << std::setprecision(8);
    std::cout << "Computing HNF decomposition over " << config.grid_length_x << "x" << config.grid_length_y << " grid." << std::endl;
    hnf::full_grid_induced_decomposition(
        config.decomposer, istream, ostream, 
        config.show_indecomp_statistics, 
        config.show_runtime_statistics, 
        config.dynamic_grid, 
        config.is_decomposed, 
        config.grid_length_x, 
        config.grid_length_y, 
        config.grassmann_value
    );
    
    return true;
}

void output_base_change_statistics(const ProgramConfig& config) {
    if (!config.decomposer.config.save_base_change) return;
    
    int total_row_ops = 0;
    for (auto& base_change : config.decomposer.base_changes) {
        total_row_ops += base_change->performed_row_ops.size();
    }
    
    if (config.decomposer.config.show_info) {
        std::cout << "Basechange: Performed " << total_row_ops << " row operations in total." << std::endl;
    }
}

void write_output(const std::ostringstream& ostream, const FileInfo& file_info, const std::string& output_string) {
    std::string output_file_path;
    hnf::write_to_file(ostream, output_file_path, file_info.input_directory, 
        file_info.file_without_extension, file_info.extension, output_string);
}

int main(int argc, char** argv) {
    ProgramConfig config;
    initialize_decomposer_config(config.decomposer);
    
    if (argc < 2) {
        if (!handle_interactive_input(argc, argv)) {
            return 1;
        }
    }
    
    if (!parse_command_line(argc, argv, config)) {
        return 0; // Help/version shown or error occurred
    }
    
    FileInfo file_info = resolve_input_file(argc, argv, config.test_files, config.is_decomposed);
    
    if (!config.test_files) {
        std::ostringstream ostream;
        if (!process_input_file(file_info, config, ostream)) {
            return 1;
        }
        
        output_base_change_statistics(config);
        
        if (config.write_output) {
            write_output(ostream, file_info, config.output_string);
        }
    }
    
    return 0;
}