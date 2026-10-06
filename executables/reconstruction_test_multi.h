#ifndef RECONSTRUCTION_TEST_MULTI_H
#define RECONSTRUCTION_TEST_MULTI_H

#include "PHY/types.h"
#include "PHY/TOOLS/tools_defs.h"
#include <omp.h> // [필수] OpenMP 헤더 추가

// Max safety sizes
#define MAX_SYM_LEN 4096 

// [수정] i7-13700F는 16코어지만 안전하게 32로 설정
#define MAX_THREADS 32 

typedef struct {
    int is_initialized;
    
    // --- 1. Signal Buffers (Thread-Specific) ---
    // [수정] MAX_ANTENNAS 대신 MAX_THREADS 사용
    double* real_raw[MAX_THREADS];
    double* imag_raw[MAX_THREADS];
    double* y_r[MAX_THREADS];
    double* y_i[MAX_THREADS];
    double* rec_r[MAX_THREADS];
    double* rec_i[MAX_THREADS];
    
    // --- Scratchpads (Thread-Specific) ---
    double* ds[MAX_THREADS];
    double* rec[MAX_THREADS];
    double* res[MAX_THREADS];
    double* bl[MAX_THREADS];
    double* egN[MAX_THREADS];
    double* s1[MAX_THREADS];
    double* s2[MAX_THREADS];
    double* tmp[MAX_THREADS];

} recon_context_t;

extern recon_context_t g_ctx;

void init_reconstruction_fast(void);
void run_fast_real_experiment(c16_t* buffer, int num_samples, int antenna_idx);
//void apply_rotation(RU_t* ru, int start_idx, int length, uint64_t timestamp_rx, int sym_offset, double f_offset);
//void apply_derotation(RU_t* ru, int start_idx, int length, uint64_t timestamp_rx, int sym_offset, double f_offset);
#endif
