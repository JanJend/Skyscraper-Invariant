#include "hnf.hpp"
#include <cassert>

int main() {
    using namespace hnf;
    // The unit square interval: one generator, two boundary relations,
    // and their one syzygy in degree (1,1). Its area is exactly one.
    R2Mat P(2, 1, {{0}, {0}}, {{1,0}, {0,1}}, {{0,0}});
    R2PModule M(P);
    M.compute_projective_resolution();
    const auto before_d1 = M.presentation();
    const auto before_d2 = M.projective_resolution()[1];
    auto working = module_presentation(M);
    working.minimize_variant();
    assert(M.projective_resolution().size() == 2);
    assert(M.presentation().data == before_d1.data);
    assert(M.presentation().col_degrees == before_d1.col_degrees);
    assert(M.projective_resolution()[1].data == before_d2.data);
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

    // A relation at (1,1) leaves unbounded strips. Both grid algorithms
    // must retain their support cutoff while preserving the supplied module.
    R2PModule strip(R2Mat(1, 1, {{0}}, {{1,1}}, {{0,0}}));
    strip.compute_projective_resolution();
    for (bool smart : {false, true}) {
        vec<R2PModule> modules{strip};
        vec<R2Mat> matrices{strip.presentation()};
        aida::AIDA_functor from_module, from_matrix;
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
        assert(modules[0].has_complete_projective_resolution());
    }
}
