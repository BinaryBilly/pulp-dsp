/* =====================================================================
 * Project:      PULP DSP Library
 * Title:        plp_mat_scale_stride_f32p_xpulpv2.c
 * Description:  parallel 32-bit floating-point strided matrix scale for XPULPV2
 *
 * $Date:        6. July 2020
 * $Revision:    V0
 *
 * Target Processor: PULP cores
 * ===================================================================== */
/*
 * Copyright (C) 2020 ETH Zurich and University of Bologna.
 *
 * Author: Tibor Schneider, ETH Zurich
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
 */

#include "plp_math.h"

/**
  @ingroup MatScaleStride
 */

/**
  @addtogroup MatScaleStrideKernels
  @{
 */

/**
   @brief Parallel strided matrix scale of 32-bit floating-point matrices kernel for XPULPV2
   extension.
   @param[in]  args      pointer to plp_mat_scale_stride_instance_f32 struct initialized by
   plp_mat_scale_stride_f32_parallel
   @return     none
*/

void plp_mat_scale_stride_f32p_xpulpv2(void *args) {

    int core_id = rt_core_id();

    plp_mat_scale_stride_instance_f32 *a = (plp_mat_scale_stride_instance_f32 *)args;

    const float *__restrict__ pSrc = a->pSrc;
    uint32_t M = a->M;
    uint32_t N = a->N;
    uint32_t strideSrc = a->strideSrc;
    uint32_t strideDst = a->strideDst;
    float scaleFactor = a->scaleFactor;
    uint32_t nPE = a->nPE;
    float *__restrict__ pDst = a->pDst;

#define BASIC_VERSION // if used don't forget to also use the undefine at end of file
#ifdef BASIC_VERSION

    for (int m = core_id; m < M; m += nPE) {
        for (int n = 0; n < N; n++) {
            float val = pSrc[m * strideSrc + n] * scaleFactor;
            pDst[m * strideDst + n] = val;
        }
    }

#else

    // TODO: Hackathon

#endif
#undef BASIC_VERSION
}

/**
   @} end of MatScaleStrideKernels group
*/
