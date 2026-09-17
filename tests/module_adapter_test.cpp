#include "uni_b1.hpp"
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
}
