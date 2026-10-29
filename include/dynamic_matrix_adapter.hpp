/**
 * @file dynamic_matrix_adapter.hpp
 * @brief Coordinate and grid matrices at Skyscraper's two-dimensional geometry boundary.
 */
#pragma once

#include <grlina/graded_linalg.hpp>
#include <grlina/matrices/coord.hpp>
#include <grlina/matrices/grid.hpp>
#include <grlina/storage/csc.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>

namespace hnf {

using R2Mat = graded_linalg::GridGradedSparseMatrix<double, 2, int>;
using R2PModule = graded_linalg::Module<R2Mat>;
using R2CSCMat = graded_linalg::R2GradedSparseMatrix<int, graded_linalg::CSCStorage<int>>;
using R2CSCModule = graded_linalg::Module<R2CSCMat>;

template <class T>
struct is_adapted_matrix : std::false_type {};

template <class Scalar, std::size_t Dimension, class Index, class Derived>
struct is_adapted_matrix<graded_linalg::CoordGradedSparseMatrix<Scalar, Dimension, Index, Derived>>
    : std::true_type {};

template <class Scalar, std::size_t Dimension, class Index>
struct is_adapted_matrix<graded_linalg::GridGradedSparseMatrix<Scalar, Dimension, Index>>
    : std::true_type {};

template <class Index, class Storage>
struct is_adapted_matrix<graded_linalg::R2GradedSparseMatrix<Index, Storage>>
    : std::true_type {};

template <class T>
struct is_adapted_module : std::false_type {};

template <class Matrix>
struct is_adapted_module<graded_linalg::Module<Matrix>> : is_adapted_matrix<Matrix> {};
template <class T>
inline constexpr bool needs_editable_presentation_v =
    (is_adapted_matrix<std::decay_t<T>>::value && !std::is_same_v<std::decay_t<T>, R2Mat>) ||
    is_adapted_module<std::decay_t<T>>::value;

template <class Matrix>
inline void require_two_parameters(const Matrix& matrix) {
    if (std::tuple_size<typename Matrix::degree_type>::value != 2)
        throw std::invalid_argument("The skyscraper invariant requires two grading parameters");
}

namespace detail {
template <class Number>
int skyscraper_index(Number value) {
    if constexpr (std::is_signed_v<Number>)
        if (value < 0)
            throw std::length_error("Skyscraper requires nonnegative matrix indices");
    if (static_cast<std::uintmax_t>(value) >
        static_cast<std::uintmax_t>(std::numeric_limits<int>::max()))
        throw std::length_error(
            "Skyscraper's int core requires ranks and grid sizes at most INT_MAX");
    return static_cast<int>(value);
}

template <class Scalar>
double skyscraper_coordinate(Scalar value) {
    const double converted = static_cast<double>(value);
    bool exact = std::isfinite(converted);
    if constexpr (std::is_integral_v<Scalar>) {
        if constexpr (std::numeric_limits<Scalar>::digits > std::numeric_limits<double>::digits) {
            // An integer is exactly representable when its odd significand fits.
            Scalar odd = value;
            while (odd != 0 && odd % 2 == 0)
                odd /= 2;
            constexpr Scalar limit = Scalar{1} << std::numeric_limits<double>::digits;
            exact = exact && odd < limit;
            if constexpr (std::is_signed_v<Scalar>)
                exact = exact && odd > -limit;
        }
    } else {
        exact = exact && static_cast<Scalar>(converted) == value;
    }
    if (!exact)
        throw std::invalid_argument(
            "Skyscraper's geometry requires coordinates exactly representable as double");
    return converted;
}
} // namespace detail

// ponytail: retain the existing double/int geometry core; template its geometry
// if coordinates outside double precision must be supported.
template <class Matrix, std::enable_if_t<is_adapted_matrix<Matrix>::value, int> = 0>
R2Mat editable_presentation(const Matrix& matrix) {
    require_two_parameters(matrix);
    if constexpr (std::is_same_v<Matrix, R2Mat>) {
        return matrix;
    } else {
        R2Mat result(detail::skyscraper_index(matrix.get_num_cols()),
                     detail::skyscraper_index(matrix.get_num_rows()));
        graded_linalg::array<int> columns;
        columns.reserve(matrix.data.size());
        for (const auto& column : matrix.data) {
            graded_linalg::vec<int> entries;
            entries.reserve(column.size());
            for (auto row : column)
                entries.push_back(detail::skyscraper_index(row));
            columns.push_back(std::move(entries));
        }
        result.assign_data(std::move(columns));
        if constexpr (graded_linalg::matrix_grid_backed_v<Matrix>) {
            // Keep every axis value: unused points can still define support
            // bounds or give meaning to degrees in later kernel computations.
            for (std::size_t axis = 0; axis < 2; ++axis) {
                detail::skyscraper_index(matrix.grid_[axis].size());
                result.grid_[axis].clear();
                for (auto value : matrix.grid_[axis])
                    result.grid_[axis].push_back(detail::skyscraper_coordinate(value));
            }
            auto degrees = [](const auto& source) {
                graded_linalg::vec<graded_linalg::base_degree<R2Mat>> converted;
                converted.reserve(source.size());
                for (auto degree : source)
                    converted.push_back(
                        {detail::skyscraper_index(degree[0]), detail::skyscraper_index(degree[1])});
                return converted;
            };
            result.col_degrees = degrees(matrix.col_degrees);
            result.row_degrees = degrees(matrix.row_degrees);
        } else {
            // Direct-coordinate matrices have no axes to preserve. Construct
            // an exact grid embedding from their physical row and column grades.
            auto degrees = [](const auto& source) {
                graded_linalg::vec<R2Mat::GridD> converted;
                converted.reserve(source.size());
                for (auto degree : source)
                    converted.push_back({detail::skyscraper_coordinate(degree[0]),
                                         detail::skyscraper_coordinate(degree[1])});
                return converted;
            };
            result.set_real_degrees(degrees(matrix.col_degrees), degrees(matrix.row_degrees));
        }
        result.compatibly_sorted = matrix.compatibly_sorted;
        return result;
    }
}

template <class Matrix, std::enable_if_t<is_adapted_matrix<Matrix>::value, int> = 0>
R2Mat editable_presentation(const graded_linalg::Module<Matrix>& module) {
    return editable_presentation(module.presentation());
}

// Convert each source once at the entry boundary. The grid algorithms keep
// their support bounds and scratch state in these editable presentations.
template <class Container>
graded_linalg::vec<R2Mat> editable_presentations(const Container& sources) {
    graded_linalg::vec<R2Mat> working;
    working.reserve(sources.size());
    for (const auto& source : sources)
        working.push_back(editable_presentation(source));
    return working;
}

template <class Degree>
graded_linalg::r2degree geometry_degree(const Degree& degree) {
    if (degree.size() != 2)
        throw std::invalid_argument("Skyscraper geometry requires a two-coordinate degree");
    return {degree[0], degree[1]};
}

inline graded_linalg::r2degree real_row_degree(const R2Mat& matrix, int row) {
    return geometry_degree(matrix.real_row_degree(row));
}

inline graded_linalg::r2degree real_col_degree(const R2Mat& matrix, int column) {
    return geometry_degree(matrix.real_col_degree(column));
}

inline graded_linalg::pair<graded_linalg::r2degree> geometry_bounds(const R2Mat& matrix) {
    require_two_parameters(matrix);
    const auto bounds = matrix.bounding_box();
    return {geometry_degree(bounds.first), geometry_degree(bounds.second)};
}

// Subspaces use the parent's embedding: the stored indices have meaning only
// together with its grids, including unused coordinates retained by kernels.
inline R2Mat graded_subspace(const R2Mat& parent, const graded_linalg::SparseMatrix<int>& basis) {
    require_two_parameters(parent);
    auto result = parent.empty_like(basis.get_num_cols(), basis.get_num_rows());
    result.assign_data(basis.data);
    result.row_degrees = parent.row_degrees;
    result.col_degrees = graded_linalg::vec<graded_linalg::base_degree<R2Mat>>(
        basis.get_num_cols(), parent.row_degrees[0]);
    result.compatibly_sorted = parent.compatibly_sorted;
    return result;
}

} // namespace hnf
