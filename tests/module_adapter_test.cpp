#define GRLINA_CSC_COMPILE_WARNINGS 0
#define GRLINA_CSC_RUNTIME_WARNINGS 0
#include "hnf.hpp"
#include <cassert>
#include <cstdint>

namespace {
using namespace hnf;

template<class Matrix>
Matrix coordinate_input(const R2Mat& reference) {
    using Degree = typename Matrix::degree_type;
    using Scalar = typename Degree::value_type;
    using Index = graded_linalg::matrix_index_t<Matrix>;
    Matrix result(reference.get_num_cols(), reference.get_num_rows());
    graded_linalg::array<Index> columns;
    for (const auto& column : reference.data) columns.emplace_back(column.begin(), column.end());
    result.assign_data(std::move(columns));
    auto degrees = [](const auto& source) {
        graded_linalg::vec<Degree> converted;
        for (const auto& degree : source) converted.push_back({static_cast<Scalar>(degree[0]), static_cast<Scalar>(degree[1])});
        return converted;
    };
    graded_linalg::GradedMatrixIO<Matrix>::assign_degrees(result,
        degrees(reference.real_col_degrees()), degrees(reference.real_row_degrees()));
    if constexpr (graded_linalg::matrix_grid_backed_v<Matrix>) {
        typename Matrix::Grid grids;
        for (std::size_t axis = 0; axis < 2; ++axis)
            grids[axis].assign(reference.grid_[axis].begin(), reference.grid_[axis].end());
        result.reindex_grid(std::move(grids));
    }
    return result;
}

void same_factors(const vec<HN_factors>& actual, const vec<HN_factors>& expected) {
    assert(actual.size() == expected.size());
    for (std::size_t i = 0; i < actual.size(); ++i) {
        assert(actual[i].size() == expected[i].size());
        for (std::size_t j = 0; j < actual[i].size(); ++j) {
            assert(actual[i][j].slope_value == expected[i][j].slope_value);
            assert(actual[i][j].d1.data == expected[i][j].d1.data);
            assert(actual[i][j].d1.real_col_degrees() == expected[i][j].d1.real_col_degrees());
            assert(actual[i][j].d1.real_row_degrees() == expected[i][j].d1.real_row_degrees());
        }
    }
}

template<class Matrix>
void accepts_coordinate_type() {
    R2Mat reference(4, 2, {{0},{0},{1},{1}},
        {{2.5,0},{0,4},{1.25,0},{0,4}}, {{0,0},{0,0}});
    if constexpr (graded_linalg::matrix_grid_backed_v<Matrix>) reference.include_real_degrees({{8,9}});
    const Matrix input = coordinate_input<Matrix>(reference);
    const auto before = input;
    auto subspaces = sparse_seperated_grassmannians<int>(2);
    vec<HN_factors> expected, actual;
    skyscraper_invariant(reference, expected, subspaces, {{0,0},{5,5}});
    skyscraper_invariant(input, actual, subspaces, {{0,0},{5,5}});
    same_factors(actual, expected);
    assert(Uni_B1(input).area() == Uni_B1(reference).area());

    R2PModule resolution(reference);
    resolution.compute_projective_resolution();
    vec<Matrix> maps;
    for (const auto& map : resolution.proj_res_.differentials()) maps.push_back(coordinate_input<Matrix>(map));
    const auto module = Module<Matrix>::from_projective_resolution(ChainComplex<Matrix>(maps), ResolutionCompleteness::complete);
    actual.clear();
    skyscraper_invariant(module, actual, subspaces, {{0,0},{5,5}});
    same_factors(actual, expected);
    Uni_B1 resolved(module);
    assert(resolved.d2.data == resolution.proj_res_[1].data);
    assert(resolved.d2.real_col_degrees() == resolution.proj_res_[1].real_col_degrees());
    const vec<Matrix> matrices{input};
    const vec<Module<Matrix>> modules{module};
    same_factors(skyscraper_invariant_sum(matrices, subspaces, {{0,0},{5,5}}), expected);
    actual.clear();
    skyscraper_invariant_sum_append(modules, actual, subspaces, {{0,0},{5,5}});
    same_factors(actual, expected);
    for (bool smart : {false,true}) {
        vec<R2Mat> working{reference};
        aida::AIDA_functor<> baseline, from_matrix, from_module;
        baseline.config.show_info = baseline.config.progress = false;
        from_matrix.config.show_info = from_matrix.config.progress = false;
        from_module.config.show_info = from_module.config.progress = false;
        std::ostringstream expected_output, matrix_output, module_output;
        if (smart) {
            process_summands_smart_grid(baseline, expected_output, 3, 3, working);
            process_summands_smart_grid(from_matrix, matrix_output, 3, 3, matrices);
            process_summands_smart_grid(from_module, module_output, 3, 3, modules);
        } else {
            process_summands_fixed_grid(baseline, expected_output, 3, 3, working);
            process_summands_fixed_grid(from_matrix, matrix_output, 3, 3, matrices);
            process_summands_fixed_grid(from_module, module_output, 3, 3, modules);
        }
        assert(matrix_output.str() == expected_output.str());
        assert(module_output.str() == expected_output.str());
    }
    assert(input.data == before.data && input.col_degrees == before.col_degrees && input.row_degrees == before.row_degrees);
    if constexpr (graded_linalg::matrix_grid_backed_v<Matrix>) assert(input.grid_ == before.grid_);
    for (std::size_t i = 0; i < maps.size(); ++i) {
        const auto& map = module.proj_res_[i];
        assert(map.data == maps[i].data && map.col_degrees == maps[i].col_degrees && map.row_degrees == maps[i].row_degrees);
        if constexpr (graded_linalg::matrix_grid_backed_v<Matrix>) assert(map.grid_ == maps[i].grid_);
    }
}
} // namespace

int main() {
    using namespace hnf;
    accepts_coordinate_type<CoordGradedSparseMatrix<double,2,int>>();
    accepts_coordinate_type<R2CSCMat>();
    accepts_coordinate_type<GridGradedSparseMatrix<float,2,std::int64_t>>();
    // Scalar == Index still uses a real grid, including non-index coordinates.
    using IntegerGrid = GridGradedSparseMatrix<std::int64_t,2,std::int64_t>;
    IntegerGrid integer_interval(2,1,{{0},{0}},{{1,0},{0,1}},{{0,0}});
    integer_interval.include_real_degrees({{-2,3}});
    assert(Uni_B1(integer_interval).area() == 1.0);
    const auto large = std::int64_t{1} << 60;
    integer_interval.include_real_degrees({{large,3}});
    assert(editable_presentation(integer_interval).grid(0).back() == static_cast<double>(large));
    integer_interval.include_real_degrees({{large + 1,3}});
    bool narrowing_refused = false;
    try { Uni_B1 unsupported(integer_interval); }
    catch (const std::invalid_argument&) { narrowing_refused = true; }
    assert(narrowing_refused); // Preserve distinct coordinates instead of merging them.
    if constexpr (std::numeric_limits<long double>::digits > std::numeric_limits<double>::digits) {
        CoordGradedSparseMatrix<long double,2,int> precise(0,1);
        precise.row_degrees[0] = {1.0L + std::numeric_limits<long double>::epsilon(),0};
        narrowing_refused = false;
        try { Uni_B1 unsupported(precise); }
        catch (const std::invalid_argument&) { narrowing_refused = true; }
        assert(narrowing_refused);
    }
    // The unit square interval: one generator, two boundary relations,
    // and their one syzygy in degree (1,1). Its area is exactly one.
    R2Mat P(2, 1, {{0}, {0}}, {{1,0}, {0,1}}, {{0,0}});
    R2PModule M(P);
    M.compute_projective_resolution();
    const auto before_d1 = M.presentation();
    const auto before_d2 = M.proj_res_[1];
    auto working = module_presentation(M);
    working.sort_compatibly();
    working.minimize_variant();
    assert(M.proj_res_.size() == 2);
    assert(M.presentation().data == before_d1.data);
    assert(M.presentation().col_degrees == before_d1.col_degrees);
    assert(M.proj_res_[1].data == before_d2.data);
    Uni_B1 interval(M);
    assert(interval.d2.data == before_d2.data);
    assert(interval.area() == 1.0);
    Uni_B1 shortcut(M.presentation(), true);
    assert(shortcut.d2.get_num_cols() == 1);
    assert(shortcut.d2.data == graded_linalg::array<int>({{0,1}}));
    assert(shortcut.area() == 1.0);
    R2Mat free(0, 1, {}, {}, {{0,0}});
    Uni_B1 free_interval(free, true);
    assert(free_interval.d2.get_num_cols() == 0);

    // Instantiate AIDA's array-degree interface, which its CLI does not use.
    aida::multipers_interface_input<int> input{{{1,0}, {0,1}}, {{0,0}}, P.data};
    aida::AIDA_functor<> decomposer;
    decomposer.config.show_info = decomposer.config.progress = false;
    auto output = decomposer.multipers_interface(input);
    assert(output.summands.size() == 1);
    auto expected_columns = input.col_degrees;
    auto actual_columns = output.summands[0].col_degrees;
    std::sort(expected_columns.begin(), expected_columns.end());
    std::sort(actual_columns.begin(), actual_columns.end());
    assert(actual_columns == expected_columns);
    assert(output.summands[0].row_degrees == input.row_degrees);
    assert(output.summands[0].matrix == P.data);

    // Real geometry must use the embedding, not the integer degree positions.
    R2PModule rectangle(R2Mat(2, 1, {{0}, {0}}, {{6.5,3}, {2.5,7}}, {{2.5,3}}));
    rectangle.compute_projective_resolution();
    assert(rectangle.proj_res_.size() == 2);
    assert((rectangle.proj_res_[0] * rectangle.proj_res_[1]).is_zero());
    Uni_B1 rectangle_factor(rectangle);
    assert(rectangle_factor.area() == 16.0);
    assert(real_col_degree(rectangle_factor.d2, 0) == r2degree({6.5,7}));
    auto last_map = rectangle.proj_res_[1];
    const auto next_kernel = last_map.graded_kernel();
    assert(next_kernel.get_num_cols() == 0);
    assert(next_kernel.grid_ == last_map.grid_);

    // The two independent rectangles have areas ten and five. Finding their
    // filtration exercises submodule kernels and a quotient between factors.
    R2Mat sum(4, 2, {{0}, {0}, {1}, {1}},
              {{2.5,0}, {0,4}, {1.25,0}, {0,4}}, {{0,0}, {0,0}});
    auto subspaces = sparse_seperated_grassmannians<int>(2);
    vec<HN_factors> factors;
    skyscraper_invariant(sum, factors, subspaces, {{0,0}, {5,5}});
    assert(factors.size() == 1 && factors[0].size() == 2);
    assert(factors[0][0].slope_value == 5.0);
    assert(factors[0][1].slope_value == 2.5);

    // CSC inputs cross the editing boundary once; physical degrees and every
    // stored resolution map remain unchanged in the caller's objects.
    R2Mat with_unused_grid = sum;
    with_unused_grid.include_real_degrees({{8,9}});
    const auto compressed = coordinate_input<R2CSCMat>(with_unused_grid);
    const auto compressed_before = compressed;
    vec<HN_factors> vector_factors, compressed_factors;
    skyscraper_invariant(with_unused_grid, vector_factors, subspaces, {{0,0}, {5,5}});
    skyscraper_invariant(compressed, compressed_factors, subspaces, {{0,0}, {5,5}});
    same_factors(compressed_factors, vector_factors);
    R2PModule resolved_vector(with_unused_grid);
    resolved_vector.compute_projective_resolution();
    vec<R2CSCMat> compressed_maps;
    for (const auto& map : resolved_vector.proj_res_.differentials())
        compressed_maps.push_back(coordinate_input<R2CSCMat>(map));
    const R2CSCModule compressed_module = R2CSCModule::from_projective_resolution(
        ChainComplex<R2CSCMat>(compressed_maps), ResolutionCompleteness::complete);
    compressed_factors.clear();
    skyscraper_invariant(compressed_module, compressed_factors, subspaces, {{0,0}, {5,5}});
    same_factors(compressed_factors, vector_factors);
    Uni_B1 compressed_resolution(compressed_module);
    assert(compressed_resolution.d2.data == resolved_vector.proj_res_[1].data);
    assert(compressed_resolution.d2.real_col_degrees() == compressed_maps[1].col_degrees);
    assert(Uni_B1(compressed).area() == Uni_B1(with_unused_grid).area());
    for (bool smart : {false, true}) {
        // R2 CSC stores physical degrees, so its grid bounds come from used degrees.
        vec<R2Mat> editable{sum};
        const vec<R2CSCMat> matrices{compressed};
        const vec<R2CSCModule> modules{compressed_module};
        aida::AIDA_functor<> vector_decomposer, matrix_decomposer, module_decomposer;
        vector_decomposer.config.show_info = vector_decomposer.config.progress = false;
        matrix_decomposer.config.show_info = matrix_decomposer.config.progress = false;
        module_decomposer.config.show_info = module_decomposer.config.progress = false;
        std::ostringstream vector_output, matrix_output, module_output;
        if (smart) {
            process_summands_smart_grid(vector_decomposer, vector_output, 3, 3, editable);
            process_summands_smart_grid(matrix_decomposer, matrix_output, 3, 3, matrices);
            process_summands_smart_grid(module_decomposer, module_output, 3, 3, modules);
        } else {
            process_summands_fixed_grid(vector_decomposer, vector_output, 3, 3, editable);
            process_summands_fixed_grid(matrix_decomposer, matrix_output, 3, 3, matrices);
            process_summands_fixed_grid(module_decomposer, module_output, 3, 3, modules);
        }
        assert(matrix_output.str() == vector_output.str());
        assert(module_output.str() == vector_output.str());
    }
    assert(compressed.data == compressed_before.data);
    assert(compressed.col_degrees == compressed_before.col_degrees);
    assert(compressed.row_degrees == compressed_before.row_degrees);
    for (std::size_t i = 0; i < compressed_maps.size(); ++i) {
        const auto& map = compressed_module.proj_res_[i];
        assert(map.data == compressed_maps[i].data);
        assert(map.col_degrees == compressed_maps[i].col_degrees);
        assert(map.row_degrees == compressed_maps[i].row_degrees);
    }

    // Baseline from the original matrix implementation on stable_example_paper:
    // one two-dimensional factor, with slope 2.42 in this padded square.
    R2Mat stable_example(5, 2, {{0}, {1}, {0,1}, {1}, {0}},
                         {{0,2}, {0,3}, {1,1}, {2,0}, {3,0}}, {{0,0}, {0,0}});
    factors.clear();
    skyscraper_invariant(stable_example, factors, subspaces, {{0,0}, {3.3,3.3}});
    assert(factors.size() == 1 && factors[0].size() == 1);
    assert(factors[0][0].d1.get_num_rows() == 2);
    assert(std::abs(factors[0][0].slope_value - 2.42) < 1e-12);

    bool refused = false;
    try {
        Uni_B1 unsupported(GridGradedSparseMatrix<double, 3, int>(0, 0));
    } catch (const std::invalid_argument&) {
        refused = true;
    }
    assert(refused);

    // Exercise the AIDA boundary inside both grid algorithms with a genuine
    // two-generator presentation; output must retain the fractional endpoint.
    for (bool smart : {false, true}) {
        vec<R2Mat> matrices{sum};
        aida::AIDA_functor<> grid_decomposer;
        grid_decomposer.config.show_info = grid_decomposer.config.progress = false;
        std::ostringstream grid_output;
        if (smart)
            process_summands_smart_grid(grid_decomposer, grid_output, 3, 3, matrices);
        else
            process_summands_fixed_grid(grid_decomposer, grid_output, 3, 3, matrices);
        assert(grid_output.str().find("(1.25;") != std::string::npos);
    }

    // A relation at (1,1) leaves unbounded strips. Both grid algorithms
    // must retain their support cutoff while preserving the supplied module.
    R2PModule strip(R2Mat(1, 1, {{0}}, {{1,1}}, {{0,0}}));
    strip.compute_projective_resolution();
    for (bool smart : {false, true}) {
        vec<R2PModule> modules{strip};
        vec<R2Mat> matrices{strip.presentation()};
        aida::AIDA_functor<> from_module, from_matrix;
        from_module.config.show_info = from_matrix.config.show_info = false;
        from_module.config.progress = from_matrix.config.progress = false;
        std::ostringstream module_output, matrix_output;
        if (smart) {
            process_summands_smart_grid(from_module, module_output, 3, 3, modules);
            process_summands_smart_grid(from_matrix, matrix_output, 3, 3, matrices);
        } else {
            process_summands_fixed_grid(from_module, module_output, 3, 3, modules);
            process_summands_fixed_grid(from_matrix, matrix_output, 3, 3, matrices);
        }
        assert(module_output.str() == matrix_output.str());
        assert(modules[0].presentation().data == strip.presentation().data);
        assert(modules[0].presentation().col_degrees == strip.presentation().col_degrees);
        assert(modules[0].res_complete_);
    }
}
