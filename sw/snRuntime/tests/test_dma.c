#include "printf.h"
#include "omp.h"
#include <stdarg.h>
#include "perf_cnt.h"
#include <snrt.h>

#define NI 32
#define NJ 32
#define NK 64
#define TILE 8
#define I_TILING ((NI%TILE==0) ? NI/TILE : NI/TILE+1) 
#define J_TILING ((NJ%TILE==0) ? NJ/TILE : NJ/TILE+1)
#define K_TILING ((NK%(TILE*2)==0) ? NK/(TILE*2) : NK/(TILE*2)+1)

const double *A_loc[2]={0x50000000, 0x50004000};
const double *B_loc[2]={0x50008000, 0x5000c000};
double *C_loc[2]={0x50010000, 0x50014000};
const uint32_t A_phys_=0xe0000360;
const uint32_t B_phys_=0xe0088360;
const uint32_t C_phys_=0xe0010360;
double alpha_=0.0000023;
double beta_=0.0000193;


void gemm_4xVL(double *c, const double *a, const double *b,
                 const unsigned int m_start, const unsigned int m_end,
                 const unsigned int N, const unsigned int P,
                 const double alpha, const double beta, const uint32_t first_iter);
int loop();

int main(){

    for(int i=0;i<1;i++){
        loop();
    }

    return 0;
}

int loop(){

    uint32_t volatile *clock=0x7800c000;

    //profiling
    uint32_t total_time = 0;
    // uint32_t total_inner_dma_wait_time = 0;
    // uint32_t total_outer_dma_wait_time = 0;
    uint32_t total_inner_dma_req_time = 0;
    uint32_t total_inner_dma_req = 0;
    uint32_t start_parall = read_csr(mcycle);


    uint32_t cluster_idx = snrt_cluster_idx();
    uint32_t global_core_idx = snrt_global_core_idx();
    uint32_t cluster_core_idx = snrt_cluster_core_idx();
    
    //repetitions over k
    uint32_t k_iters = K_TILING;
    uint32_t dim_k_stride = NK / k_iters;

    //repetitions over i
    uint32_t i_iters = I_TILING;
    uint32_t dim_i_stride = NI / i_iters;
    uint32_t start_index_I = dim_i_stride/2 * cluster_idx; // index of each cluster in the I dimension
    uint32_t core_size_I = (dim_i_stride/4);
    uint32_t loc_start_index_i = (cluster_core_idx==0) ? 0 : (core_size_I);

    //repetitions over j
    uint32_t j_iters = J_TILING;
    uint32_t dim_j_stride = NJ / j_iters;

    //double buffering
    uint32_t curr_A_B =0;
    uint32_t curr_C =0;

    //values of C dma copy
    uint32_t size_C_transfer = dim_j_stride*sizeof(double);

    //values of A dma copy
    uint32_t size_A_transfer = NK*sizeof(double)/k_iters;

    //values of B dma copy
    uint32_t size_B_transfer = dim_j_stride*sizeof(double);

    //first copy of the double buffering
    if(cluster_core_idx==0){
        snrt_dma_start_2d_wideptr(  (uint64_t)C_loc[0], //dest
                                (uint64_t)C_phys_ + start_index_I*NJ*sizeof(double), //src
                                size_C_transfer, //size
                                size_C_transfer, //dst_stride
                                NJ*sizeof(double), //src_stride
                                (dim_i_stride/2) //repeat
        );

        snrt_dma_start_2d_wideptr(  (uint64_t)A_loc[0], //dest
                                    (uint64_t)A_phys_ + start_index_I*NK*sizeof(double), //src
                                    size_A_transfer, //size
                                    size_A_transfer, //dst_stride
                                    NK*sizeof(double), //src_stride
                                    (dim_i_stride/2) //repeat
        );

        snrt_dma_start_2d_wideptr(  (uint64_t)B_loc[0], //dest
                                    (uint64_t)B_phys_, //src
                                    size_B_transfer, //size
                                    size_B_transfer, //dst_stride
                                    NJ*sizeof(double), //src_stride
                                    dim_k_stride //repeat
        );   
        // uint32_t start_outer_dma_wait = read_csr(mcycle);   
        snrt_dma_wait_all();   
        // uint32_t end_outer_dma_wait = read_csr(mcycle);
        // total_outer_dma_wait_time += end_outer_dma_wait-start_outer_dma_wait;

    }
    snrt_cluster_hw_barrier();

    for(int i=0;i<i_iters;i++){
        for(int j=0;j<j_iters;j++){

            double *C_curr[1]; 
            double *C_next[1]; 
            if(curr_C == 0){
                C_curr[0] = C_loc[0];
                C_next[0] = C_loc[1];
                curr_C = 1;
            }else{
                C_curr[0] = C_loc[1];
                C_next[0] = C_loc[0];
                curr_C = 0;
            }

            
            if(cluster_core_idx==0){ 
                // uint32_t start_inner_dma_req = read_csr(mcycle);
                if(j!=j_iters-1){
                    snrt_dma_start_2d_wideptr(  (uint64_t)C_next[0], //dest
                                        (uint64_t)C_phys_ + (start_index_I*NJ + i*dim_i_stride*NJ + (j+1)*dim_j_stride)*sizeof(double), //src
                                        size_C_transfer, //size
                                        size_C_transfer, //dst_stride
                                        NJ*sizeof(double), //src_stride
                                        (dim_i_stride/2) //repeat
                    );
                    // total_inner_dma_req++;

                }else if(i!=i_iters-1){
                    snrt_dma_start_2d_wideptr(  (uint64_t)C_next[0], //dest
                                            (uint64_t)C_phys_ + (start_index_I*NJ + (i+1)*dim_i_stride*NJ)*sizeof(double), //src
                                            size_C_transfer, //size
                                            size_C_transfer, //dst_stride
                                            NJ*sizeof(double), //src_stride
                                            (dim_i_stride/2) //repeat
                    );
                    // total_inner_dma_req++;
                }
                // uint32_t end_inner_dma_req = read_csr(mcycle);
                // total_inner_dma_req_time += end_inner_dma_req-start_inner_dma_req;

            }

            snrt_cluster_hw_barrier();  


            for(int k=0;k<k_iters;k++){
        
                double *A_curr[1];  
                double *A_next[1];
                double *B_curr[1]; 
                double *B_next[1];

                if(curr_A_B == 0){
                    A_curr[0] = A_loc[0];
                    A_next[0] = A_loc[1];
                    B_curr[0] = B_loc[0];
                    B_next[0] = B_loc[1];
                    curr_A_B = 1;
                }else{
                    A_curr[0] = A_loc[1];
                    A_next[0] = A_loc[0];
                    B_curr[0] = B_loc[1];
                    B_next[0] = B_loc[0];
                    curr_A_B = 0;
                }

                if(cluster_core_idx==0){
                    if(k!=0 || k_iters==1){
                        // uint32_t start_dma_wait = read_csr(mcycle);
                        snrt_dma_wait_all();
                        // uint32_t end_dma_wait = read_csr(mcycle);
                        // total_inner_dma_wait_time += end_dma_wait-start_dma_wait;
                    }
                    uint32_t start_inner_dma_req = read_csr(mcycle);
                    if(k!=k_iters-1){
                        snrt_dma_start_2d_wideptr(  (uint64_t)A_next[0], //dest
                                                    (uint64_t)A_phys_ + (start_index_I*NK + (k+1)*dim_k_stride + i*dim_i_stride*NK)*sizeof(double), //src
                                                    size_A_transfer, //size
                                                    size_A_transfer, //dst_stride
                                                    NK*sizeof(double), //src_stride
                                                    (dim_i_stride/2) //repeat
                        );
                        total_inner_dma_req++;

                    }else if(j!=j_iters-1){
                        snrt_dma_start_2d_wideptr(  (uint64_t)A_next[0], //dest
                                                    (uint64_t)A_phys_ + (start_index_I*NK + i*dim_i_stride*NK)*sizeof(double), //src
                                                    size_A_transfer, //size
                                                    size_A_transfer, //dst_stride
                                                    NK*sizeof(double), //src_stride
                                                    (dim_i_stride/2) //repeat
                        );
                        total_inner_dma_req++;
                    }else if(i!=i_iters-1){
                        snrt_dma_start_2d_wideptr(  (uint64_t)A_next[0], //dest
                                                    (uint64_t)A_phys_ + (start_index_I*NK + (i+1)*dim_i_stride*NK)*sizeof(double), //src
                                                    size_A_transfer, //size
                                                    size_A_transfer, //dst_stride
                                                    NK*sizeof(double), //src_stride
                                                    (dim_i_stride/2) //repeat
                        );
                        total_inner_dma_req++;
                    }

                    if(k!=k_iters-1){
                    
                        snrt_dma_start_2d_wideptr(  (uint64_t)B_next[0], //dest
                                                    (uint64_t)B_phys_ + ((k+1)*dim_k_stride*NJ + j*dim_j_stride)*sizeof(double), //src
                                                    size_B_transfer, //size
                                                    size_B_transfer, //dst_stride
                                                    NJ*sizeof(double), //src_stride
                                                    dim_k_stride //repeat
                        );    
                        total_inner_dma_req++;

                    }else if(j!=j_iters-1){

                        snrt_dma_start_2d_wideptr(  (uint64_t)B_next[0], //dest
                                                    (uint64_t)B_phys_ +  (j+1)*dim_j_stride*sizeof(double), //src
                                                    size_B_transfer, //size
                                                    size_B_transfer, //dst_stride
                                                    NJ*sizeof(double), //src_stride
                                                    dim_k_stride //repeat
                        );  
                        total_inner_dma_req++;

                    }else if(i!=i_iters-1){
                        snrt_dma_start_2d_wideptr(  (uint64_t)B_next[0], //dest
                                                    (uint64_t)B_phys_, //src
                                                    size_B_transfer, //size
                                                    size_B_transfer, //dst_stride
                                                    NJ*sizeof(double), //src_stride
                                                    dim_k_stride //repeat
                        );  
                        total_inner_dma_req++;

                    }
                    
                    uint32_t end_inner_dma_req = read_csr(mcycle);
                    total_inner_dma_req_time += end_inner_dma_req-start_inner_dma_req;
                }

                snrt_cluster_hw_barrier();

                uint32_t first_iter = (k==0) ? 1 : 0;
                uint32_t start = read_csr(mcycle);

                gemm_4xVL(C_curr[0], A_curr[0], B_curr[0], 
                                    loc_start_index_i, loc_start_index_i+core_size_I, dim_k_stride, dim_j_stride, alpha_, beta_,first_iter);
                

                snrt_cluster_hw_barrier();
                uint32_t end = read_csr(mcycle);
                total_time += end-start;

            }
            if(cluster_core_idx==0){
                snrt_dma_start_2d_wideptr( (uint64_t)C_phys_ + (start_index_I*NJ + i*dim_i_stride*NJ + j*dim_j_stride)*sizeof(double), //dest
                                            (uint64_t)C_curr[0], //src
                                            size_C_transfer, //size
                                            NJ*sizeof(double), //dst_stride
                                            size_C_transfer, //src_stride
                                            (dim_i_stride/2) //repeat

                );
            }
        }
    }

    if(cluster_core_idx==0){
        // uint32_t start_last_dma_wait = read_csr(mcycle);
        snrt_dma_wait_all();  
        // uint32_t end_last_dma_wait = read_csr(mcycle);
        // total_outer_dma_wait_time += end_last_dma_wait-start_last_dma_wait;
    }

    uint32_t end_parall = read_csr(mcycle);

    if(snrt_global_core_idx()==0)*clock=total_inner_dma_req_time;
    return 0;
}

void gemm_4xVL(double *c, const double *a, const double *b,
                 const unsigned int m_start, const unsigned int m_end,
                 const unsigned int N, const unsigned int P,
                 const double alpha, const double beta, const uint32_t first_iter) {
                    

    unsigned int p = 0;
    const unsigned int P_striding = P*sizeof(double);
    const unsigned int N_striding = N*sizeof(double);
    while (p < P) {
        // Calculate the vl
        size_t gvl;
        asm volatile("vsetvli %[gvl], %[vl], e64, m4, ta, ma"
                    : [gvl] "=r"(gvl)
                    : [vl] "r"(P - p));

        const double *b_ = b + p;
        double *c_ = c + p;

        for (unsigned int m = m_start; m < m_end; m += 4) {
            const double *a_ = a + m * N;
            const double *a__ = a_;

            asm volatile("vle64.v v16, (%0);" ::"r"(b_));
            const double *b__;
            // b__ = b_ + P;
            asm volatile("add %0, %1, %2" : "+r"( b__) :"r"( b_), "r"(P_striding));

            double *c__ = c_ + m * P;


            double t0, t1, t2, t3;


            asm volatile(   "fld   %[t0], (%[a__])              \n"         // t0 = *a__;
                            "add   %[a__], %[a__], %[N_striding]\n"         // a__ += N;
                            "fld    %[t1], (%[a__])             \n"
                            "add   %[a__], %[a__], %[N_striding]\n"
                            "fld    %[t2], (%[a__])             \n"
                            "add   %[a__], %[a__], %[N_striding]\n"
                            "fld    %[t3], (%[a__])             \n"
                    :   [t0] "+f"(t0), [t1] "+f"(t1),[t2] "+f"(t2), [t3] "+f"(t3), [a__] "+r"(a__)
                    : [N_striding]"r"(N_striding)
            );


            unsigned int n = 0;
            // uint32_t start = read_csr(mcycle);
            while (n < N_striding) {
                // a__ = a_ + ++n;
                asm volatile("addi %[n], %[n], %[incr]\n"
                            "add %[a__], %[a_], %[n]\n"
                    : [n] "+r"(n), [a__] "+r"(a__)
                    : [incr]"i"(sizeof(double)), [a_]"r"(a_)
                );


                asm volatile("vle64.v v20, (%0);" ::"r"(b__));
                // b__ += P;
                asm volatile("add %0, %0, %1" : "+r"( b__) : "r"(P_striding));


                if (n == sizeof(double)) {

                    asm volatile(   "vfmul.vf v0, v16, %[t0]            \n"
                                    "fld    %[t0], (%[a__])             \n"
                                    "add   %[a__], %[a__], %[N_striding]\n"
                                    "vfmul.vf v4, v16, %[t1]            \n"
                                    "fld    %[t1], (%[a__])             \n"
                                    "add   %[a__], %[a__], %[N_striding]\n"
                                    "vfmul.vf v8, v16, %[t2]            \n"
                                    "fld    %[t2], (%[a__])             \n"
                                    "add   %[a__], %[a__], %[N_striding]\n"
                                    "vfmul.vf v12, v16, %[t3]           \n"
                                    "fld    %[t3], (%[a__])             \n"

                        : [t0] "+f"(t0), [t1] "+f"(t1),[t2] "+f"(t2), [t3] "+f"(t3), [a__] "+r"(a__)
                        : [N_striding] "r"(N_striding)
                    );
                


                } else {

                    asm volatile(   "vfmacc.vf v0, %[t0], v16           \n"
                                    "fld    %[t0], (%[a__])             \n"
                                    "add   %[a__], %[a__], %[N_striding]\n"
                                    "vfmacc.vf v4, %[t1], v16           \n"
                                    "fld    %[t1], (%[a__])             \n"
                                    "add   %[a__], %[a__], %[N_striding]\n"
                                    "vfmacc.vf v8, %[t2], v16           \n"
                                    "fld    %[t2], (%[a__])             \n"
                                    "add   %[a__], %[a__], %[N_striding]\n"
                                    "vfmacc.vf v12, %[t3], v16          \n"
                                    "fld    %[t3], (%[a__])             \n"

                        : [t0] "+f"(t0), [t1] "+f"(t1),[t2] "+f"(t2), [t3] "+f"(t3), [a__] "+r"(a__)
                        : [N_striding] "r"(N_striding)
                    );
                }

                // a__ = a_ + ++n;
                asm volatile("addi %[n], %[n], %[incr]\n"
                            "add %[a__], %[a_], %[n]\n"
                    : [n] "+r"(n), [a__] "+r"(a__)
                    : [incr]"i"(sizeof(double)), [a_]"r"(a_)
                );

                if (n == N_striding)
                    break;

                asm volatile("vle64.v v16, (%0);" ::"r"(b__));
                // b__ += P;
                asm volatile("add %0, %0, %1" : "+r"( b__) : "r"(P_striding));

                asm volatile(   "vfmacc.vf v0, %[t0], v20           \n"
                                "fld    %[t0], (%[a__])             \n"
                                "add   %[a__], %[a__], %[N_striding]\n"
                                "vfmacc.vf v4, %[t1], v20           \n"
                                "fld    %[t1], (%[a__])             \n"
                                "add   %[a__], %[a__], %[N_striding]\n"
                                "vfmacc.vf v8, %[t2], v20           \n"
                                "fld    %[t2], (%[a__])             \n"
                                "add   %[a__], %[a__], %[N_striding]\n"
                                "vfmacc.vf v12, %[t3], v20          \n"
                                "fld    %[t3], (%[a__])             \n"

                    : [t0] "+f"(t0), [t1] "+f"(t1),[t2] "+f"(t2), [t3] "+f"(t3), [a__] "+r"(a__)
                    : [N_striding] "r"(N_striding)
                );
            }
            // uint32_t end = read_csr(mcycle);
            // if(snrt_global_core_idx()==0){
            //     snrt_printf("\n%x,N %u, Inner cycle: %u\n\r",N,  end-start);
            // }

            //last accumulation
            asm volatile("vfmacc.vf v0, %0, v20" ::"f"(t0));
            asm volatile("vfmacc.vf v4, %0, v20" ::"f"(t1));
            asm volatile("vfmacc.vf v8, %0, v20" ::"f"(t2));
            asm volatile("vfmacc.vf v12, %0, v20" ::"f"(t3));

            double *c___ = c__;
            
            //multiply by alpha and hide load of c vectors inside
            asm volatile(   "vle64.v v16, (%[c])            \n" 
                            "add %[c], %[c], %[P]           \n"
                            "vfmul.vf v0, v0, %[alpha]      \n"
                            "vle64.v v20, (%[c])            \n"
                            "add %[c], %[c], %[P]           \n"
                            "vfmul.vf v4, v4, %[alpha]      \n"
                            "vle64.v v24, (%[c])            \n"
                            "add %[c], %[c], %[P]           \n"
                            "vfmul.vf v8, v8, %[alpha]      \n"
                            "vle64.v v28, (%[c])            \n"
                            "vfmul.vf v12, v12, %[alpha]    \n"

            :[c]"+r"(c___)
            :[alpha]"f"(alpha), [P]"r"(P_striding) );


            if(first_iter==1){

                asm volatile(   "vfmul.vf v16, v16, %[beta] \n"
                                "vfmul.vf v20, v20, %[beta] \n"
                                "vfmul.vf v24, v24, %[beta] \n"
                                "vfmul.vf v28, v28, %[beta] \n"

                                "vfadd.vv v0, v0, v16       \n"
                                "vse64.v v0, (%[c])         \n"
                                "add %[c], %[c], %[P]       \n"
                                
                                "vfadd.vv v4, v4, v20       \n"
                                "vse64.v v4, (%[c])         \n"
                                "add %[c], %[c], %[P]       \n"
                                
                                "vfadd.vv v8, v8, v24       \n"
                                "vse64.v v8, (%[c])         \n"
                                "add %[c], %[c], %[P]       \n"

                                "vfadd.vv v12, v12, v28     \n"
                                "vse64.v v12, (%[c])        \n"
                : [c]"+r"(c__)
                :[beta] "f"(beta), [P]"r"(P_striding)
                );

            }else{
                //add to c and store result
                asm volatile(   
                                "vfadd.vv v0, v0, v16       \n"
                                "vse64.v v0, (%[c])         \n"
                                "add %[c], %[c], %[P]       \n"
                                
                                "vfadd.vv v4, v4, v20       \n"
                                "vse64.v v4, (%[c])         \n"
                                "add %[c], %[c], %[P]       \n"
                                
                                "vfadd.vv v8, v8, v24       \n"
                                "vse64.v v8, (%[c])         \n"
                                "add %[c], %[c], %[P]       \n"

                                "vfadd.vv v12, v12, v28     \n"
                                "vse64.v v12, (%[c])        \n"
                : [c]"+r"(c__)
                : [P]"r"(P_striding)
                );                
            }
        
        }

        p += gvl;
    }
}