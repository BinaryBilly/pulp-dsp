/* =====================================================================
 * Project:      PULP DSP Library
 * Title:        plp_cmplx_conj_f32_xpulpv2.c
 * Description:  16-bit integer complex dot product glue code
 *
 * $Date:        29. June 2020
 * $Revision:    V0
 *
 * Target Processor: PULP cores
 * ===================================================================== */
/*
 * Copyright (C) 2019 ETH Zurich and Ubiversity of Bologna.
 *
 * Author: Hanna Mueller, ETH Zurich
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * Licensed under the Apache License, Version 2.0 (the License); you may
 * not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS, WITHOUT
 * WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Notice: project inspired by ARM CMSIS DSP and parts of source code
  ported and adopted for RISC-V PULP platform from ARM CMSIS DSP
 released under Copyright (C) 2010-2019 ARM Limited or its affiliates
  with Apache-2.0.
 */

#include "plp_math.h"

/**
  @ingroup groupCmplxMath
 */

/**
  @defgroup cmplx_dot_prod Complex Dot Product
  Computes the dot product of two complex vectors.
  The vectors are multiplied element-by-element and then summed.
  The <code>pSrcA</code> points to the first complex input vector and
  <code>pSrcB</code> points to the second complex input vector.
  <code>numSamples</code> specifies the number of complex samples
  and the data in each array is stored in an interleaved fashion
  (real, imag, real, imag, ...).
  Each array has a total of <code>2*numSamples</code> values.
  The underlying algorithm is used:
  <pre>
  realResult = 0;
  imagResult = 0;
  for (n = 0; n < numSamples; n++) {
      realResult += pSrcA[(2*n)+0] * pSrcB[(2*n)+0] - pSrcA[(2*n)+1] * pSrcB[(2*n)+1];
      imagResult += pSrcA[(2*n)+0] * pSrcB[(2*n)+1] + pSrcA[(2*n)+1] * pSrcB[(2*n)+0];
  }
  </pre>
  There are separate functions for floating point, integer, and fixed point 32- 16- 8-bit data
  types.
 */

/**
  @addtogroup cmplx_dot_prod
  @{
 */

/**
  @brief         16-bit integer complex dot product.
  @param[in]     pSrcA       points to the first input vector
  @param[in]     pSrcB       points to the second input vector
  @param[in]     numSamples  number of samples in each vector
  @param[out]    realResult  real part of the result returned here
  @param[out]    imagResult  imaginary part of the result returned here
  @return        none
 */

void plp_cmplx_dot_prod_i16_xpulpv2(const int16_t *pSrcA,
                                    const int16_t *pSrcB,
                                    uint32_t numSamples,
                                    int16_t *realResult,
                                    int16_t *imagResult) {
    uint32_t i;
    int32_t real_sum = 0;
    int32_t imag_sum = 0;
    uint32_t loopEnd = numSamples/2;

    for (i = 0; i < loopEnd; i++)
    {
      // int32_t r_a1 = pSrcA[4*i];
      // int32_t i_a1 = pSrcA[4*i+1];
      // int32_t r_a2 = pSrcA[4*i+2];
      // int32_t i_a2 = pSrcA[4*i+3];
      // int32_t r_b1 = pSrcB[4*i];
      // int32_t i_b1 = pSrcB[4*i+1];
      // int32_t r_b2 = pSrcB[4*i+2];
      // int32_t i_b2 = pSrcB[4*i+3];

      v2s ab = *((v2s *)&(pSrcA[4*i]));
      v2s ef = *((v2s *)&(pSrcA[4*i+2]));
      v2s cd = *((v2s *)&(pSrcB[4*i]));
      v2s gh = *((v2s *)&(pSrcB[4*i+2]));

      // real_sum += ab[0]*cd[0] - ab[1]*cd[1] + ef[0]*gh[0] - ef[1]*gh[1];
      v2s bf = __builtin_shuffle(ab,ef,(v2s){1,3});
      v2s dh = __builtin_shuffle(cd,gh,(v2s){1,3});
      v2s ae = __builtin_shuffle(ab,ef,(v2s){0,2});
      v2s cg = __builtin_shuffle(cd,gh,(v2s){0,2});
      // real_sum += -bf[0]*dh[0] - bf[1]*dh[1] + ae[0]*cg[0] + ae[1]*cg[1];
      real_sum = __SUMDOTP2(ae,cg,real_sum);
      real_sum -= __DOTP2(bf,dh);

      // imag_sum += ab[0]*cd[1] + ab[1]*cd[0] + ef[0]*gh[1] + ef[1]*gh[0];
      v2s dc = __builtin_shuffle(cd,cd,(v2s){1,0});
      v2s hg = __builtin_shuffle(gh,gh,(v2s){1,0});
      // imag_sum += ab[0]*dc[0] + ab[1]*dc[1] + ef[0]*hg[0] + ef[1]*hg[1];
      imag_sum = __SUMDOTP2(ab, dc, imag_sum);
      imag_sum = __SUMDOTP2(ef, hg, imag_sum);

    }

    if(numSamples&1){
      int32_t r_a1 = pSrcA[4*i];
      int32_t i_a1 = pSrcA[4*i+1];
      int32_t r_b1 = pSrcB[4*i];
      int32_t i_b1 = pSrcB[4*i+1];
      real_sum += r_a1*r_b1 - i_a1*i_b1;
      imag_sum += r_a1*i_b1 + i_a1*r_b1;
    }


    *realResult = real_sum;
    *imagResult = imag_sum;
}
/**
  @} end of cmplx_dot_prod group
 */