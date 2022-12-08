
#include "field.h"

#define DEC_IDX                                                  \
    int blk_i = (blockIdx.z * gridDim.y * gridDim.x) +           \
                (blockIdx.y * gridDim.x) +                       \
                blockIdx.x;                                      \
                                                                 \
    int thd_i = (blk_i * blockDim.z * blockDim.y * blockDim.x) + \
                (threadIdx.z * blockDim.y * blockDim.x) +        \
                (threadIdx.y * blockDim.x) + threadIdx.x;

typedef struct combined_t
{
    Fr table_value;
    Fr a_minus;
    ulong r_next;
    ulong r_prev;
} combined_t;
static_assert(sizeof(combined_t) == 80);

extern "C" __global__ void compute_evaluate_h_lookups_codeblock(
    Fr *values,
    const combined_t *combined_data_in,
    const Fr *product_coset_list,
    const Fr *permuted_input_coset_list,
    const Fr *permuted_table_coset_list,
    const Fr *l0,
    const Fr *l_active_row,
    const Fr *l_last,
    const Fr *y_beta_gamma_one,
    const int lookup_count,
    const int array_size)
{

    DEC_IDX

    const Fr y = y_beta_gamma_one[0];
    const Fr beta = y_beta_gamma_one[1];
    const Fr gamma = y_beta_gamma_one[2];
    const Fr one = y_beta_gamma_one[3];

    const Fr l0_ = l0[thd_i];
    const Fr l_active_row_ = l_active_row[thd_i];
    const Fr l_last_ = l_last[thd_i];

    Fr value = values[thd_i];

    for (int n = 0; n < lookup_count; n++)
    {

        const int list_offset = array_size * n;

        const Fr table_value = combined_data_in[list_offset + thd_i].table_value;
        const Fr a_minus_s = combined_data_in[list_offset + thd_i].a_minus;
        const ulong r_next = combined_data_in[list_offset + thd_i].r_next;
        const ulong r_prev = combined_data_in[list_offset + thd_i].r_prev;

        const Fr *product_coset = &product_coset_list[list_offset];
        const Fr *permuted_input_coset = &permuted_input_coset_list[list_offset];
        const Fr *permuted_table_coset = &permuted_table_coset_list[list_offset];

        // l_0(X) * (1 - z(X)) = 0
        value = value * y + ((one - product_coset[thd_i]) * l0_);
        // l_last(X) * (z(X)^2 - z(X)) = 0
        value = value * y +
                ((product_coset[thd_i] *
                      product_coset[thd_i] -
                  product_coset[thd_i]) *
                 l_last_);
        // (1 - (l_last(X) + l_blind(X))) * (
        //   z(\omega X) (a'(X) + \beta) (s'(X) + \gamma)
        //   - z(X) (\theta^{m-1} a_0(X) + ... + a_{m-1}(X) + \beta) (\theta^{m-1}
        //     s_0(X) + ... + s_{m-1}(X) + \gamma)
        // ) = 0
        value = value * y +
                ((product_coset[r_next] *
                      (permuted_input_coset[thd_i] + beta) *
                      (permuted_table_coset[thd_i] + gamma) -
                  product_coset[thd_i] * table_value) *
                 l_active_row_);
        // Check that the first values in the permuted input expression and permuted
        // fixed expression are the same.
        // l_0(X) * (a'(X) - s'(X)) = 0
        value = value * y + (a_minus_s * l0_);
        // Check that each value in the permuted lookup input expression is either
        // equal to the value above it, or the value at the same index in the
        // permuted table expression.
        // (1 - (l_last + l_blind)) * (a′(X) − s′(X))⋅(a′(X) − a′(\omega^{-1} X)) =
        // 0
        value = value * y +
                (a_minus_s *
                 (permuted_input_coset[thd_i] -
                  permuted_input_coset[r_prev]) *
                 l_active_row_);
    }

    values[thd_i] = value;
}