

#include "file_reader.hpp"
#include <iostream>
#include <fstream>
#include <sstream>  
#include <regex>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cassert>

namespace hnf {

    double safe_stod(const std::string& str, const std::string& context) {
    if (str.empty()) {
        throw std::runtime_error("Empty string in " + context);
    }
    try {
        return std::stod(str);
    } catch (const std::exception& e) {
        std::cerr << "stod failed at " << context << ": '" << str << "'" << std::endl;
        throw;
    }
}

GridData bars_from_sky(const std::string& filename) {
    std::ifstream file(filename);
    static char buffer[1 << 20]; // 1 MB buffer
    file.rdbuf()->pubsetbuf(buffer, sizeof(buffer));
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    GridData result;
    std::string line;
    
    // Line 1: Must be "HNF"
    std::getline(file, line);
    if (line.find("HNF") == std::string::npos) {
        throw std::runtime_error("First line must be 'HNF'");
    }
    
    // Line 2: Grid dimensions
    std::getline(file, line);
    std::sscanf(line.c_str(), "%d,%d", &result.n_x, &result.n_y);
    
    std::cout << "Grid dimensions for filtered landscape: " << result.n_x << " x " << result.n_y << std::endl;

    // Line 3: Lattice info - extract coordinates
    std::getline(file, line);
    std::regex coord_regex(R"(\((-?[0-9.eE+-]+),\s*(-?[0-9.eE+-]+)\))");
    auto coords_begin = std::sregex_iterator(line.begin(), line.end(), coord_regex);
    auto coords_end = std::sregex_iterator();
    
    std::vector<std::pair<double, double>> coords;
    for (auto it = coords_begin; it != coords_end; ++it) {
        double x = safe_stod((*it)[1], "lattice vector x-coordinate");
        double y = safe_stod((*it)[2], "lattice vector y-coordinate");
        coords.push_back({x, y});
    }
    
    if (coords.size() < 3) {
        throw std::runtime_error("Expected 3 coordinate pairs");
    }
    
    // Extract lattice vectors
    result.start_x = coords[0].first;
    result.start_y = coords[0].second;
    result.end_x = coords[1].first;
    result.end_y = coords[1].second;
    result.step_x = coords[2].first;
    result.step_y = coords[2].second;
    result.slope = result.step_y / result.step_x;
    
    // Initialize bars grid
    result.bars.resize(result.n_x, std::vector<std::vector<Bar>>(result.n_y));
    
    // Process grid points and stable modules
    std::pair<double, double> current_position;
    int i = -1, j = -1;
    
    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);
        
        if (line.empty()) continue;
        
        if (line.substr(0, 2) == "G,") {
            size_t first_comma = line.find(',', 2);
            size_t second_comma = line.find(',', first_comma + 1);
            size_t paren_open = line.find('(', second_comma);
            size_t coord_comma = line.find(',', paren_open);
            size_t paren_close = line.find(')', coord_comma);
            
            i = std::stoi(line.substr(2, first_comma - 2));
            j = std::stoi(line.substr(first_comma + 1, second_comma - first_comma - 1));
            current_position.first = safe_stod(line.substr(paren_open + 1, coord_comma - paren_open - 1), "grid point x-coordinate");
            current_position.second = safe_stod(line.substr(coord_comma + 1, paren_close - coord_comma - 1), "grid point y-coordinate");
        } else {

            std::istringstream iss(line);
            std::string token;
            
            // First token is slope
            std::getline(iss, token, ',');
            double theta = safe_stod(token, "slope");
            
            // Parse relations
            std::vector<std::pair<double, double>> relations;
            while (std::getline(iss, token, ',')) {
                size_t paren_open = token.find('(');
                size_t semicolon = token.find(';');
                size_t paren_close = token.find(')');
                
                if (paren_open != std::string::npos && 
                    semicolon != std::string::npos && 
                    paren_close != std::string::npos) {
                    double x = safe_stod(token.substr(paren_open + 1, semicolon - paren_open - 1), "relation x-coordinate");
                    double y = safe_stod(token.substr(semicolon + 1, paren_close - semicolon - 1), "relation y-coordinate");
                    relations.push_back({x, y});
                } else {
                    throw std::runtime_error("Malformed relation: " + token);
                }
            }
            
            if (relations.empty()) continue;
            
            // Convert to relative coordinates and find intersections
            std::vector<std::pair<double, double>> candidates;
            for (const auto& [x, y] : relations) {
                double rel_x = x - current_position.first;
                double rel_y = y - current_position.second;
                
                double int_x, int_y;
                if (rel_y <= result.slope * rel_x) {
                    int_x = rel_x;
                    int_y = result.slope * rel_x;
                } else {
                    int_y = rel_y;
                    int_x = rel_y / result.slope;
                }
                candidates.push_back({int_x, int_y});
            }
            
            // Choose minimum candidate
            auto min_candidate = *std::min_element(candidates.begin(), candidates.end());
            
            // Calculate length (magnitude)
            double length = min_candidate.first;
            result.bars[i][j].push_back({theta, length});
        }
    }
    // result.n_y -= 30; // Adjust for extra grid points
    std::cout << "Loaded landscape grid of size " << result.n_x << " x " << result.n_y << std::endl;
    return result;
}

} // namespace hnf