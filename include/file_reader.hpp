#pragma once

#ifndef FILE_READER_HPP
#define FILE_READER_HPP

#include <vector>
#include <string>

namespace hnf {
struct Bar {
    double theta;
    double length;
};

struct GridData {
    int n_x, n_y;
    double start_x, start_y, end_x, end_y, step_x, step_y;
    double slope;
    std::vector<std::vector<std::vector<Bar>>> bars;
};

// Function declarations
GridData bars_from_sky(const std::string& filename);
} // namespace hnf

#endif // FILE_READER_HPP