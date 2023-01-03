#include <arbor/fvm_types.hpp>

#include <arbor/gpu/gpu_common.hpp>
#include "matrix_common.hpp"

namespace arb {
namespace gpu {

namespace kernels {
/// GPU implementatin of Hines matrix assembly
/// Flat layout
/// For a given time step size dt
///     - use the precomputed alpha and alpha_d values to construct the diagonal
///       and off diagonal of the symmetric Hines matrix.
///     - compute the RHS of the linear system to solve
template <typename T, typename I>
__global__
void assemble_matrix_flat(
        T* d,
        T* rhs,
        const T* invariant_d,
        const T* voltage,
        const T* current,
        const T* conductivity,
        const T* cv_capacitance,
        const T* cv_area,
        const I* cv_to_cell,
        const T* dt_intdom,
        const I* cell_to_intdom,
        unsigned n)
{
    const unsigned tid = threadIdx.x + blockDim.x*blockIdx.x;

    if (tid<n) {
        auto cid = cv_to_cell[tid];
        auto dt = dt_intdom[cell_to_intdom[cid]];

        // Note: dt==0 case is expected only at the end of a mindelay/2
        // integration period, and consequently divergence is unlikely
        // to be a peformance problem.

        if (dt>0) {
            // The 1e-3 is a constant of proportionality required to ensure that the
            // conductance (gi) values have units μS (micro-Siemens).
            // See the model documentation in docs/model for more information.
            T oodt_factor = 1e-3/dt; // [1/μs]
            T area_factor = 1e-3*cv_area[tid]; // [1e-9·m²]

            auto gi = oodt_factor * cv_capacitance[tid] + area_factor*conductivity[tid]; // [μS]
            d[tid] = gi + invariant_d[tid];
            rhs[tid] = gi*voltage[tid] - area_factor*current[tid];
        }
        else {
            d[tid] = 0;
            rhs[tid] = voltage[tid];
        }
    }
}

} // namespace kernels

ARB_ARBOR_API void assemble_matrix_flat(
        arb_value_type* d,
        arb_value_type* rhs,
        const arb_value_type* invariant_d,
        const arb_value_type* voltage,
        const arb_value_type* current,
        const arb_value_type* conductivity,
        const arb_value_type* cv_capacitance,
        const arb_value_type* area,
        const arb_index_type* cv_to_cell,
        const arb_value_type* dt_intdom,
        const arb_index_type* cell_to_intdom,
        unsigned n)
{
    constexpr unsigned block_dim = 128;
    const unsigned grid_dim = impl::block_count(n, block_dim);

    kernels::assemble_matrix_flat
        <arb_value_type, arb_index_type>
        <<<grid_dim, block_dim>>>
        (d, rhs, invariant_d, voltage, current, conductivity, cv_capacitance,
         area, cv_to_cell, dt_intdom, cell_to_intdom, n);
}

} // namespace gpu
} // namespace arb
