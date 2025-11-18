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


int get_diagonal_index(int i, int j, const GridData& data);

void compute_landscape(const GridData& data, const std::string& output_filename, 
                      double theta, int k);

} // namespace hnf

#endif // LANDSCAPE_HPP