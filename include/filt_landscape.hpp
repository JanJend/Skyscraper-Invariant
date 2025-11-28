#ifndef LANDSCAPE_HPP
#define LANDSCAPE_HPP

#include <vector>
#include <string>
#include "file_reader.hpp"

namespace hnf {

struct GridPoint {
    int i, j;
    double x, y;
};

struct RGB {
    unsigned char r, g, b;
};

RGB heatmap_color(double value);

void write_landscape_png(const std::vector<std::vector<std::vector<double>>>& landscape, 
                         const std::string& output_filename);

void write_landscape(const std::vector<std::vector<double>>& landscape, const GridData& data, const std::string& output_filename,
                      const double& theta_min, const int& k);

int get_diagonal_index(int i, int j, const GridData& data);

std::vector<std::vector<std::vector<double>>> compute_landscape(const GridData& data, 
                      const double& theta, const int& k);

} // namespace hnf

#endif // LANDSCAPE_HPP