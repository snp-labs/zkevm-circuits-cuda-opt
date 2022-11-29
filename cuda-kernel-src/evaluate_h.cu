
#include "field.h"

#define DEC_IDX                                                  \
    int blk_i = (blockIdx.z * gridDim.y * gridDim.x) +           \
                (blockIdx.y * gridDim.x) +                       \
                blockIdx.x;                                      \
                                                                 \
    int thd_i = (blk_i * blockDim.z * blockDim.y * blockDim.x) + \
                (threadIdx.z * blockDim.y * blockDim.x) +        \
                (threadIdx.y * blockDim.x) + threadIdx.x;

extern "C" __global__ void compute_evaluate_h_lookups_codeblock(
    Fr *values,
    const Fr *table_values,
    const ulong *r_next_list,
    const ulong *r_prev_list,
    const Fr *a_minus_s_list,
    const Fr *product_coset,
    const Fr *permuted_input_coset,
    const Fr *permuted_table_coset,
    const Fr *l0,
    const Fr *l_active_row,
    const Fr *l_last,
    const Fr *y_beta_gamma_one)
{

    DEC_IDX

    const Fr table_value = table_values[thd_i];
    const ulong r_next = r_next_list[thd_i];
    const ulong r_prev = r_prev_list[thd_i];
    const Fr a_minus_s = a_minus_s_list[thd_i];

    const Fr y = y_beta_gamma_one[0];
    const Fr beta = y_beta_gamma_one[1];
    const Fr gamma = y_beta_gamma_one[2];
    const Fr one = y_beta_gamma_one[3];

    Fr value = values[thd_i];

    // l_0(X) * (1 - z(X)) = 0
    value = value * y + ((one - product_coset[thd_i]) * l0[thd_i]);
    // l_last(X) * (z(X)^2 - z(X)) = 0
    value = value * y + ((product_coset[thd_i] * product_coset[thd_i] - product_coset[thd_i]) * l_last[thd_i]);
    // (1 - (l_last(X) + l_blind(X))) * (
    //   z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
    //   - z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta) (\theta^{m-1}
    //     s_0(X) + ... + s_{m-1}(X) + \gamma)
    // ) = 0
    value = value * y + ((product_coset[r_next] * (permuted_input_coset[thd_i] + beta) * (permuted_table_coset[thd_i] + gamma) - product_coset[thd_i] * table_value) * l_active_row[thd_i]);
    // Check that the first values in the permuted input expression and permuted
    // fixed expression are the same.
    // l_0(X) * (a'(X) - s'(X)) = 0
    value = value * y + (a_minus_s * l0[thd_i]);
    // Check that each value in the permuted lookup input expression is either
    // equal to the value above it, or the value at the same index in the
    // permuted table expression.
    // (1 - (l_last + l_blind)) * (a′(X) − s′(X))⋅(a′(X) − a′(\omega^{-1} X)) =
    // 0
    value = value * y + (a_minus_s * (permuted_input_coset[thd_i] - permuted_input_coset[r_prev]) * l_active_row[thd_i]);

    values[thd_i] = value;
}