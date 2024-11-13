#include "printf.h"
#include "omp.h"
#include <stdarg.h>
#include "perf_cnt.h"
#include <snrt.h>
#define N_BUFFER 512 //64 N per buffer
#define K_BUFFER 8



int main(){
    uint32_t volatile *clock=0x7800c000;
    const double *a=0x50000000;
    const double *b=0x50008000;
    const double *y=0x50010000;
    const double *x=0x50012000;
    int index_start = 4*snrt_global_core_idx();
    double alpha = 0.2;
    double beta = 0.01;
    unsigned int N_stride = N_BUFFER;
    int first_iter = 0;
    uint32_t begin=0, end=0;
    uint32_t result;
    for(uint32_t iters=0;iters<2;iters++){
        begin = read_csr(mcycle);
        const double *a_ = a + index_start*N_stride;
        const double *b_ = b + index_start*N_stride;
        const double *x_ = x;
        unsigned int vl;
        unsigned int avl = N_stride;
        unsigned int N_striding = N_stride*sizeof(double);

        asm volatile("vsetvli %0, %1, e64, m2, ta, ma" : "=r"(vl) : "r"(avl));

        const double *a__ = a_;
        const double *b__ = b_;

        asm volatile( 
                        "vle64.v v0, (%[a__])\n"
                        "vle64.v v4, (%[b__])\n"
                        "add  %[a__], %[a__], %[N_striding]\n"
                        "add  %[b__], %[b__], %[N_striding]\n"
                        "vle64.v v8, (%[x_])\n"
            :[a__] "+r"(a__), [b__] "+r"(b__) 
            :[x_]"r"(x_), [N_striding] "r"(N_striding)
            : "memory"
        );

        int inner_first_iter = 1;
        do{
            if(inner_first_iter && first_iter){
            asm volatile( 
                            "vle64.v v2,  (%[a__])\n"
                            "vle64.v v6,  (%[b__])\n"
                            "sll  t0,  %[vl], 3\n"
                            "vfmul.vv v24, v0, v8\n"
                            "add  %[a_], %[a_], t0\n"
                            "add  %[b_], %[b_], t0\n"
                            "vfmul.vv v28, v4, v8\n"
                            "add  %[x_], %[x_], t0\n"
                :[a_]"+r"(a_), [b_]"+r"(b_), [x_]"+r"(x_)
                :[N_striding]"r"(N_striding), [vl]"r"(vl),[a__]"r"(a__), [b__]"r"(b__)
                :"memory", "t0"
            );
            }else{
            asm volatile( 
                            "vle64.v v2,  (%[a__])\n"
                            "vle64.v v6,  (%[b__])\n"
                            "sll  t0,  %[vl], 3\n"
                            "vfmacc.vv v24, v0, v8\n"
                            "add  %[a_], %[a_], t0\n"
                            "add  %[b_], %[b_], t0\n"
                            "vfmacc.vv v28, v4, v8\n"
                            "add  %[x_], %[x_], t0\n"
                :[a_]"+r"(a_), [b_]"+r"(b_), [x_]"+r"(x_)
                : [N_striding]"r"(N_striding), [vl]"r"(vl),[a__]"r"(a__), [b__]"r"(b__)
                :"memory", "t0"
            );
            // flops += 2*16;
            }
            avl -= vl;
            if(avl>0){
            a__= a_;
            b__= b_;
            if(inner_first_iter && first_iter){
                asm volatile( 
                            "mv %[inner_first_iter], zero\n"
                            "vle64.v v0,  (%[a__])\n"
                            "vle64.v v4, (%[b__])\n"
                            "vsetvli %[vl], %[avl], e64, m2, ta, ma\n"
                            "add  %[a__], %[a__], %[N_striding]\n"
                            "vfmul.vv v26, v2, v8\n"
                            "add  %[b__], %[b__], %[N_striding]\n"
                            "vfmul.vv v30, v6, v8\n"
                            "vle64.v v8, (%[x_])\n"
                :[vl]"=r"(vl), [a__]"+r"(a__), [b__]"+r"(b__), [x_]"+r"(x_), [inner_first_iter]"=r"(inner_first_iter)
                :[avl]"r"(avl), [N_striding]"r"(N_striding)
                : "memory"
                );
            }else{
                asm volatile( 
                            "vle64.v v0,  (%[a__])\n"
                            "vle64.v v4, (%[b__])\n"
                            "vsetvli %[vl], %[avl], e64, m2, ta, ma\n"
                            "add  %[a__], %[a__], %[N_striding]\n"
                            "vfmacc.vv v26, v2, v8\n"
                            "add  %[b__], %[b__], %[N_striding]\n"
                            "vfmacc.vv v30, v6, v8\n"
                            "vle64.v v8, (%[x_])\n"
                :[vl]"=r"(vl), [a__]"+r"(a__), [b__]"+r"(b__), [x_]"+r"(x_)
                :[avl]"r"(avl), [N_striding]"r"(N_striding)
                : "memory"
                );             
            };
            }else{
            if (inner_first_iter && first_iter) {
                asm volatile("vfmul.vv v26, v2, v8");
                asm volatile("vfmul.vv v30, v6, v8");
            } else {
                asm volatile("vfmacc.vv v26, v2, v8");
                asm volatile("vfmacc.vv v30, v6, v8");
            }
            }
        }while(avl>0);

        end = read_csr(mcycle);
        if(iters==0) result=end-begin;
        else if(end-begin<result) result = end-begin;
    }
    if(snrt_global_core_idx()==0)*clock=result;
    return 0;
}