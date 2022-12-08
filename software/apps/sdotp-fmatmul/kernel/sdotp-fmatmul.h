// Copyright 2021 ETH Zurich and University of Bologna.
//
// SPDX-License-Identifier: Apache-2.0
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

// Author: Domenic Wüthrich, ETH Zurich

#ifndef SDOTPFMATMUL_H
#define SDOTPFMATMUL_H

void matmul(__fp16 *c, const __fp16 *a, const __fp16 *b, const unsigned int M,
            const unsigned int N, const unsigned int P);

inline void matmul_2xVL(__fp16 *c, const __fp16 *a, const __fp16 *b,
                        const unsigned int m_start, const unsigned int m_end,
                        const unsigned int N, const unsigned int P,
                        const unsigned int p_start, const unsigned int p_end,
                        const unsigned int vl) __attribute__((always_inline));
inline void matmul_4xVL(__fp16 *c, const __fp16 *a, const __fp16 *b,
                        const unsigned int m_start, const unsigned int m_end,
                        const unsigned int N, const unsigned int P,
                        const unsigned int p_start, const unsigned int p_end,
                        const unsigned int vl) __attribute__((always_inline));
inline void matmul_8xVL(__fp16 *c, const __fp16 *a, const __fp16 *b,
                        const unsigned int m_start, const unsigned int m_end,
                        const unsigned int N, const unsigned int P,
                        const unsigned int p_start, const unsigned int p_end,
                        const unsigned int vl) __attribute__((always_inline));

#endif
