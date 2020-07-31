/* =====================================================================
 * Project:      PULP DSP Library
 * Title:        plp_cmplx_conj_i8_xpulpv2.c
 * Description:  32-bit float complex conjugate glue code
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
  @defgroup cmplx_conj Complex Conjugate
  Conjugates the elements of a complex data vector.
  The <code>pSrc</code> points to the source data and
  <code>pDst</code> points to the destination data where the result should be written.
  <code>numSamples</code> specifies the number of complex samples
  and the data in each array is stored in an interleaved fashion
  (real, imag, real, imag, ...).
  Each array has a total of <code>2*numSamples</code> values.
  The underlying algorithm is used:
  <pre>
  for (n = 0; n < numSamples; n++) {
      pDst[(2*n)  ] =  pSrc[(2*n)  ];    // real part
      pDst[(2*n)+1] = -pSrc[(2*n)+1];    // imag part
  }
  </pre>
  There are separate functions for floating point, integer, and fixed point 32- 16- 8-bit data
  types.
 */

/**
  @addtogroup cmplx_conj
  @{
 */

/**
  @brief         8-bit integer complex conjugate.
  @param[in]     pSrc        points to the input vector
  @param[out]    pDst        points to the output vector
  @param[in]     numSamples  number of samples in each vector
  @return        none
 */

void plp_cmplx_conj_i8_xpulpv2(const int8_t *__restrict__ pSrc,
                               int8_t *__restrict__ pDst,
                               uint32_t numSamples) {
    uint32_t i = 0;
    uint32_t loop_length = numSamples/2;
    uint32_t clean_up = numSamples & 1;
    v4s cVec;
    v4s inverter = {0,-1,0,-1};
    v4s adder = {0,1,0,1};

    for (i = 0; i < loop_length; i++)
    {
      cVec = *((v4s* )&(pSrc[4*i]));
      cVec = __EXOR4(cVec, inverter);
      cVec = __ADD4(cVec, adder);
      *((v4s *)&(pDst[4*i])) = cVec;
    }

    if(clean_up){
      pDst[i*4  ] = pSrc[i*4];
      pDst[i*4+1] = -pSrc[i*4+1];
    }

}

/**
  @} end of cmplx_conj group
 */