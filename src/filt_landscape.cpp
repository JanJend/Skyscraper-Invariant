#include "filt_landscape.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cassert>
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

namespace hnf {


RGB heatmap_color(double value) {
    // Clamp to [0,1]
    value = std::clamp(value, 0.0, 1.0);
    
    RGB color;
    if (value < 0.33) {
        // Black to red
        double t = value / 0.33;
        color.r = static_cast<unsigned char>(255 * t);
        color.g = 0;
        color.b = 0;
    } else if (value < 0.66) {
        // Red to yellow
        double t = (value - 0.33) / 0.33;
        color.r = 255;
        color.g = static_cast<unsigned char>(255 * t);
        color.b = 0;
    } else {
        // Yellow to white
        double t = (value - 0.66) / 0.34;
        color.r = 255;
        color.g = 255;
        color.b = static_cast<unsigned char>(255 * t);
    }
    return color;
}

void write_landscape_png(const std::vector<std::vector<std::vector<double>>>& landscapes,
                         const std::string& output_filename) {
    if (landscapes.empty() || landscapes[0].empty()) return;
    
    int width = landscapes[0].size();
    int height = landscapes[0][0].size();
    
    // Find min/max from k=1 landscape (index 0) for normalization
    double min_val = landscapes[0][0][0], max_val = landscapes[0][0][0];
    for (const auto& col : landscapes[0]) {
        for (double val : col) {
            min_val = std::min(min_val, val);
            max_val = std::max(max_val, val);
        }
    }
    double range = max_val - min_val;
    if (range < 1e-10) range = 1.0;
    
    // Write PNG for each k value
    for (int k_idx = 0; k_idx < landscapes.size(); k_idx++) {
        const auto& landscape = landscapes[k_idx];
        
        // Create RGB image data
        std::vector<unsigned char> image(width * height * 3);
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                double normalized = (landscape[i][j] - min_val) / range;
                RGB color = heatmap_color(normalized);
                int idx = ((height - 1 - j) * width + i) * 3;
                image[idx + 0] = color.r;
                image[idx + 1] = color.g;
                image[idx + 2] = color.b;
            }
        }
        
        // Modify filename to include k value
        std::string filename = output_filename;
        size_t ext_pos = filename.find_last_of('.');
        if (ext_pos != std::string::npos) {
            filename.insert(ext_pos, "_k" + std::to_string(k_idx + 1));
        } else {
            filename += "_k" + std::to_string(k_idx + 1);
        }
        
        stbi_write_png(filename.c_str(), width, height, 3, image.data(), width * 3);
        std::cout << "PNG landscape written to " << filename << std::endl;
    }
}

// obsolete
int get_diagonal_index(int i, int j, const GridData& data) {
    // Grid point coordinates
    double x = data.start_x + i * data.step_x;
    double y = data.start_y + j * data.step_y;
    double S = data.slope;
    
    // Line through (x, y) with slope S: Y = S*X + (y - S*x)
    double intercept = y - S * x;
    
    // Check intersection with left edge (X = start_x)
    double y_at_left = S * data.start_x + intercept;
    
    // Check intersection with bottom edge (Y = start_y)
    double x_at_bottom = (data.start_y - intercept) / S;
    
    // Grid boundaries
    double max_x = data.start_x + (data.n_x - 1) * data.step_x;
    double max_y = data.start_y + (data.n_y - 1) * data.step_y;
    
    // Determine which boundary the diagonal starts from
    bool intersects_left = (y_at_left >= data.start_y && y_at_left <= max_y);
    bool intersects_bottom = (x_at_bottom >= data.start_x && x_at_bottom <= max_x);
    
    if (intersects_bottom) {
        // Diagonal starts from bottom edge
        // Find which grid column index this corresponds to
        int col_idx = std::round((x_at_bottom - data.start_x) / data.step_x);
        return col_idx;
    } else if (intersects_left) {
        // Diagonal starts from left edge (excluding bottom-left corner)
        // Find which grid row index this corresponds to
        int row_idx = std::round((y_at_left - data.start_y) / data.step_y);
        return data.n_x + row_idx - 1;
    } else {
        throw std::runtime_error("Grid point does not lie on any valid diagonal");
    }
}

void write_landscape(const std::vector<std::vector<double>>& landscape, const GridData& data, const std::string& output_filename,
                     const double& theta_min, const int& k){
    // Write output
    std::ofstream out(output_filename);
    out << std::fixed << std::setprecision(8);
    out << "Sky Landscape " << data.n_x << " " << data.n_y << " " << k << " " << theta_min << "\n";
    for (int j = 0; j < data.n_y; j++) {
        for (int i = 0; i < data.n_x; i++) {
            out << landscape[i][j] << " ";
        }
        out << "\n";
    }
    std::cout << "Landscape written to " << output_filename << std::endl;
    out.close();
}

std::vector<std::vector<std::vector<double>>> compute_landscape(const GridData& data,
                      const double& theta_min, const int& k) {

    // Initialize landscape array
    std::vector<std::vector<std::vector<double>>> landscapes(k, std::vector<std::vector<double>>(data.n_x, std::vector<double>(data.n_y, 0.0)));
    // Process each grid point
    for (int i = 0; i < data.n_x; i++) {
        for (int j = 0; j < data.n_y; j++) {
            auto bars = data.bars[i][j];

            // Filter out bars with theta < theta_min
            bars.erase(std::remove_if(bars.begin(), bars.end(),
                    [theta_min](const Bar& bar) { return bar.theta < theta_min; }),
                    bars.end());

            if (bars.empty()) {
                continue;
            }
            
            std::sort(bars.begin(), bars.end(),
                    [](const Bar& a, const Bar& b) { return a.length > b.length; });
            
            // Process each k value
            for (int k_idx = 1; k_idx <= std::min(k, static_cast<int>(bars.size())); k_idx++) {
                const Bar& chosen = bars[k_idx - 1];
                double length = chosen.length;
                double d = length / 2.0;
                
                // Update landscape along the diagonal
                for (int t = 0; i + t < data.n_x && j + t < data.n_y; t++) {
                    double value = std::max(0.0, d - std::abs(d - t * data.step_x));
                    if (value > landscapes[k_idx - 1][i + t][j + t]) {
                        landscapes[k_idx - 1][i + t][j + t] = value;
                    }
                    if(t > 0 && value == 0.0){
                        break;
                    }
                }
            }
        }
    }
    return landscapes;
}

std::vector<std::vector<std::vector<double>>> compute_difference_landscape(const GridData& data, 
                      const double& theta,
                      const double& theta_prime, 
                     const int& k){
    assert(theta_prime <= theta);
    // Initialize landscape array
    std::vector<std::vector<std::vector<double>>> landscapes(k, std::vector<std::vector<double>>(data.n_x, std::vector<double>(data.n_y, 0.0)));
    std::vector<std::vector<std::vector<double>>> landscapes_copy(k, std::vector<std::vector<double>>(data.n_x, std::vector<double>(data.n_y, 0.0)));
    // Process each grid point
    for (int i = 0; i < data.n_x; i++) {
        for (int j = 0; j < data.n_y; j++) {
            auto bars = data.bars[i][j];
            auto bars_copy = bars;
            // Filter out bars with theta < theta_min
            bars.erase(std::remove_if(bars.begin(), bars.end(),
                    [theta](const Bar& bar) { return bar.theta < theta; }),
                    bars.end());
            // Filter out bars with theta < theta_prime
            bars_copy.erase(std::remove_if(bars_copy.begin(), bars_copy.end(),
                    [theta_prime](const Bar& bar) { return bar.theta < theta_prime; }),
                    bars_copy.end());
            if (bars.empty()) {
                continue;
            }
            
            std::sort(bars.begin(), bars.end(),
                    [](const Bar& a, const Bar& b) { return a.length > b.length; });
            std::sort(bars_copy.begin(), bars_copy.end(),
                    [](const Bar& a, const Bar& b) { return a.length > b.length; });

            // Process each k value
            for (int k_idx = 1; k_idx <= std::min(k, static_cast<int>(bars.size())); k_idx++) {
                const Bar& chosen = bars[k_idx - 1];
                double length = chosen.length;
                double d = length / 2.0;
                
                // Update landscape along the diagonal
                for (int t = 0; i + t < data.n_x && j + t < data.n_y; t++) {
                    double value = std::max(0.0, d - std::abs(d - t * data.step_x));
                    if (value > landscapes[k_idx - 1][i + t][j + t]) {
                        landscapes[k_idx - 1][i + t][j + t] = value;
                    }
                    if(t > 0 && value == 0.0){
                        break;
                    }
                }
            }
            // Now compute the landscape with theta_prime filtered bars
            for (int k_idx = 1; k_idx <= std::min(k, static_cast<int>(bars_copy.size())); k_idx++) {
                const Bar& chosen = bars_copy[k_idx - 1];
                double length = chosen.length;      
                double d = length / 2.0;
                // Update landscape along the diagonal
                for (int t = 0; i + t < data.n_x && j + t < data.n_y; t++) {
                    double value = std::max(0.0, d - std::abs(d - t * data.step_x));
                    if (value > landscapes_copy[k_idx - 1][i + t][j + t]) {
                        landscapes_copy[k_idx - 1][i + t][j + t] = value;
                    }
                    if(t > 0 && value == 0.0){
                        break;
                    }
                }
            }
        }
    }
    // Compute difference
    for(int k_idx =0; k_idx < k; k_idx++){
        for(int i =0; i < data.n_x; i++){
            for(int j =0; j < data.n_y; j++){
                assert(landscapes_copy[k_idx][i][j] >= landscapes[k_idx][i][j]);
                landscapes_copy[k_idx][i][j] -= landscapes[k_idx][i][j];
            }
        }
    }
    return landscapes_copy;
}

} // namespace hnf
