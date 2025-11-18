#include "filt_landscape.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <cassert>

namespace hnf {

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

void compute_landscape(const GridData& data, const std::string& output_filename,
                      double theta_min, int k) {

    // Initialize landscape array
    std::vector<std::vector<double>> landscape(data.n_x, std::vector<double>(data.n_y, 0.0));

    // Process each grid point
    for (int i = 0; i < data.n_x; i++) {
        for (int j = 0; j < data.n_y; j++) {
            auto bars = data.bars[i][j];

            // Filter out bars with theta < theta_min
            bars.erase(std::remove_if(bars.begin(), bars.end(),
                       [theta_min](const Bar& bar) { return bar.theta < theta_min; }),
                       bars.end());

            if (bars.size() < k) {
                continue;
            }
            std::sort(bars.begin(), bars.end(),
                      [](const Bar& a, const Bar& b) { return a.length > b.length; });
            const Bar& chosen = bars[k - 1];
            double length = chosen.length;
            double d = length / 2.0;
            
            // Update landscape along the diagonal
            for (int t = 0; i + t < data.n_x && j + t < data.n_y; t++) {
                double value = std::max(0.0, d - std::abs(d - t * data.step_x));
                if (value > landscape[i + t][j + t]) {
                    landscape[i + t][j + t] = value;
                } 
                if(t > 0 && value == 0.0){
                    break;
                }
            }
        }
    }
    
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

} // namespace hnf
