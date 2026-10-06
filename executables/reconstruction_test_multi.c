#include "reconstruction_test_multi.h"
#include "common/utils/LOG/log.h"
#include "PHY/defs_common.h"
#include "PHY/defs_nr_common.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h> // [필수] OpenMP 함수 사용

// Global Context Definition
recon_context_t g_ctx = {0};

// --- CONFIGURATION ---
#define LAMBDA          0.05      // Threshold
#define RECON_N         2            // Difference Order
static const double q15_scale = 1.0 / 32768.0;

// --- FAST HELPER FUNCTIONS (Inlined) ---
static inline double H_lambda(double input, double lambda) {
    double two_lambda = 2.0 * lambda;
    double frac_part = (input / two_lambda + 0.5) - floor(input / two_lambda + 0.5);
    return two_lambda * (frac_part - 0.5);
}

static inline double RD(double x, double lambda) {
    double two_lambda = 2.0 * lambda;
    return round(x / two_lambda) * two_lambda;
}

// --- INITIALIZATION (Called Once) ---
void init_reconstruction_fast(void) {
    if (g_ctx.is_initialized) return;
    
    // [수정] 시스템의 최대 쓰레드 개수를 확인하되, MAX_THREADS(헤더파일 정의)를 넘지 않도록 제한
    int max_sys_threads = omp_get_max_threads();
    if (max_sys_threads > MAX_THREADS) {
        LOG_W(PHY, "[RECON] Warning: System threads (%d) > MAX_THREADS (%d). Clamping.\n", max_sys_threads, MAX_THREADS);
        max_sys_threads = MAX_THREADS;
    }

    LOG_I(PHY, "[RECON] Allocating THREAD-SAFE memory for %d threads...\n", max_sys_threads);
    int pad_len = MAX_SYM_LEN + 5;

    // [중요] 안테나 개수가 아니라 '쓰레드 개수'만큼 메모리를 할당합니다.
    for (int i = 0; i < max_sys_threads; i++) {
        // Data Buffers
        g_ctx.real_raw[i] = (double*)malloc(MAX_SYM_LEN * sizeof(double));
        g_ctx.imag_raw[i] = (double*)malloc(MAX_SYM_LEN * sizeof(double));
        g_ctx.y_r[i]      = (double*)malloc(MAX_SYM_LEN * sizeof(double));
        g_ctx.y_i[i]      = (double*)malloc(MAX_SYM_LEN * sizeof(double));
        g_ctx.rec_r[i]    = (double*)malloc(MAX_SYM_LEN * sizeof(double));
        g_ctx.rec_i[i]    = (double*)malloc(MAX_SYM_LEN * sizeof(double));

        // Scratchpads (Allocated per THREAD)
        g_ctx.ds[i]  = (double*)malloc(pad_len * sizeof(double));
        g_ctx.rec[i] = (double*)malloc(pad_len * sizeof(double));
        g_ctx.res[i] = (double*)malloc(pad_len * sizeof(double));
        g_ctx.bl[i]  = (double*)malloc(pad_len * sizeof(double));
        g_ctx.egN[i] = (double*)malloc(pad_len * sizeof(double));
        g_ctx.s1[i]  = (double*)malloc(pad_len * sizeof(double));
        g_ctx.s2[i]  = (double*)malloc(pad_len * sizeof(double));
        g_ctx.tmp[i] = (double*)malloc(pad_len * sizeof(double));
    }

    g_ctx.is_initialized = 1;
}

// --- CORE RECONSTRUCTION (Optimized, No Malloc) ---
// [수정] 마지막 인자로 ant_idx 대신 tid(쓰레드 ID)를 받습니다.
static void reconstruct_fast_core(double* y, double* recovered, int size, double lambda, int N, double beta_g, int tid) {
    
    // [핵심] 안테나 번호가 아닌 '쓰레드 번호(tid)'로 메모리를 가져옵니다.
    // 그래야 동시 실행 시 충돌이 없습니다.
    double* ds  = g_ctx.ds[tid];
    double* rec = g_ctx.rec[tid];
    double* res = g_ctx.res[tid];
    double* bl  = g_ctx.bl[tid];
    double* egN = g_ctx.egN[tid];
    double* s1  = g_ctx.s1[tid];
    double* s2  = g_ctx.s2[tid];
    double* tmp = g_ctx.tmp[tid];

    // 1. Compute N-th Difference
    memcpy(ds, y, size * sizeof(double));
    int curr_len = size;
    
    for(int k=0; k<N; k++) {
        for(int i=0; i<curr_len-1; i++) tmp[i] = ds[i+1] - ds[i];
        tmp[curr_len-1] = 0;
        curr_len--;
        memcpy(ds, tmp, curr_len * sizeof(double));
    }

    // 2. Modulo & Residual
    for(int i=0; i<curr_len; i++) {
        rec[i] = H_lambda(ds[i], lambda);
        res[i] = RD(rec[i] - ds[i], lambda);
    }

    // 3. Backward Loop
    memcpy(bl, res, curr_len * sizeof(double));
    memcpy(egN, res, curr_len * sizeof(double));
    int bl_len = curr_len;
    int J = (int)ceil(6.0 * beta_g / lambda);

    for (int k = 1; k < N; k++) {
        // Cumsum 1
        s1[0] = 0; for(int i=0; i<bl_len; i++) s1[i+1] = s1[i] + egN[i];
        // Cumsum 2
        s2[0] = 0; for(int i=0; i<bl_len+1; i++) s2[i+1] = s2[i] + s1[i];

        // Kappa
        int idx_1 = 1; 
        int idx_J = (J + 1 >= bl_len + 2) ? bl_len + 1 : J + 1;
        double Kn = floor((-s2[idx_J] + s2[idx_1]) / (12.0 * beta_g) + 0.5);

        // Integrate & Correct
        tmp[0] = 0; for(int i=0; i<bl_len; i++) tmp[i+1] = tmp[i] + bl[i];
        bl_len++;
        
        double corr = 2.0 * Kn * lambda;
        for(int i=0; i<bl_len; i++) bl[i] = RD(tmp[i], lambda) + corr;
        
        memcpy(egN, bl, bl_len * sizeof(double));
    }

    // 4. Final Integration
    tmp[0] = 0; for(int i=0; i<bl_len; i++) tmp[i+1] = tmp[i] + bl[i];
    
    for(int i=0; i<size; i++) recovered[i] = tmp[i] + y[i];
}

// --- MAIN RUN FUNCTION (Called inside Parallel Region) ---
void run_fast_real_experiment(c16_t* buffer, int num_samples, int antenna_idx) {
    
    // Safety check (보통 ru_thread 루프 시작 전에 init을 호출하지만, 혹시 모르니 체크)
    if (!g_ctx.is_initialized) init_reconstruction_fast();

    // [핵심 변경] 현재 이 함수를 실행하고 있는 쓰레드의 ID를 가져옵니다.
    // antenna_idx는 로그를 찍을 때만 쓰고, 실제 메모리 접근에는 tid를 씁니다.
    int tid = omp_get_thread_num();

    // 안전장치: 쓰레드 ID가 할당된 범위를 넘어가면 실행 중지
    if (tid >= MAX_THREADS) {
        // 너무 많은 쓰레드가 생성된 경우 (설정 오류)
        return; 
    }

    // Use THREAD-LOCAL Buffers (Not Antenna-Specific!)
    double* r_raw = g_ctx.real_raw[tid];
    double* i_raw = g_ctx.imag_raw[tid];
    double* y_r   = g_ctx.y_r[tid];
    double* y_i   = g_ctx.y_i[tid];
    double* rec_r = g_ctx.rec_r[tid];
    double* rec_i = g_ctx.rec_i[tid];

    double max_amp = 0.0;
    //double sum_r = 0.0;
    //double sum_i = 0.0;
    
    // 1. Convert & Measure
    for (int k = 0; k < num_samples; k++) {
        r_raw[k] = (double)buffer[k].r * q15_scale;
        i_raw[k] = (double)buffer[k].i * q15_scale;
        
        //sum_r += r_raw[k];
        //sum_i += i_raw[k];
    }

    // DC 오프셋 계산 (전체 평균)
    //double dc_r = sum_r / num_samples;
    //double dc_i = sum_i / num_samples;
    /*
    // [Pass 2] DC 제거 및 Max Amplitude 재계산
    // 여기서 500 같은 DC를 빼주지 않으면 Reconstruction 알고리즘이 동작하지 않습니다.
    for (int k = 0; k < num_samples; k++) {
        // 핵심: 원본에서 DC 성분을 뺍니다.
        r_raw[k] -= dc_r; 
        i_raw[k] -= dc_i;
        
        // DC가 빠진 "순수 신호"를 기준으로 진폭 측정
        if(fabs(r_raw[k]) > max_amp) max_amp = fabs(r_raw[k]);
        if(fabs(i_raw[k]) > max_amp) max_amp = fabs(i_raw[k]);
    }
    */
    
    // 신호가 너무 작을 때를 대비한 안전장치
    if (max_amp < LAMBDA) max_amp = LAMBDA * 1.1;
    
    /*
    double current_max = 0.0;
    
    // 1) 현재 안테나로 들어온 신호 중 가장 큰 값(Peak) 찾기
    for (int k = 0; k < num_samples; k++) {
        if (fabs(r_raw[k]) > current_max) current_max = fabs(r_raw[k]);
        if (fabs(i_raw[k]) > current_max) current_max = fabs(i_raw[k]);
    }

    // 2) 목표 한계값 설정 (0.09)
    // 이 값을 넘어가면 Reconstruction 알고리즘이 동작하지 않는다고 판단하여 줄입니다.
    double safety_limit = 0.09; 

    // 3) 스케일링 팩터 계산
    // 신호가 0.09보다 클 때만 동작 (작으면 1.0 그대로 둠)
    double scale_factor = 1.0;
    
    if (current_max > safety_limit) {
        // 예: 신호가 0.5면 -> 0.09 / 0.5 = 0.18배로 축소
        scale_factor = safety_limit / current_max;
        
        // (선택사항) 디버깅용 로그: 얼마나 줄어들었는지 확인
        // if (antenna_idx == 1) {
        //     LOG_D(PHY, "[Limiter] Ant 1 too huge (%.4f). Scaling by %.4f to fit %.2f\n", 
        //           current_max, scale_factor, safety_limit);
        // }
    }

    // 4) 스케일링 적용 (비율 유지하며 축소)
    if (scale_factor != 1.0) {
        for (int k = 0; k < num_samples; k++) {
            r_raw[k] *= scale_factor;
            i_raw[k] *= scale_factor;
        }
        
        // 줄어든 값으로 max_amp 변수 업데이트 (매우 중요)
        // 이걸 안 하면 뒤쪽 로직에서 옛날(큰) 값을 써서 오동작할 수 있음
        max_amp = 0.0;
        for (int k = 0; k < num_samples; k++) {
             if (fabs(r_raw[k]) > max_amp) max_amp = fabs(r_raw[k]);
             if (fabs(i_raw[k]) > max_amp) max_amp = fabs(i_raw[k]);
        }
    }
    */
    // =====================================================================
    // [여기서부터는 기존 코드 유지]
    // 2. Fold ... 
    // =====================================================================
    for (int k = 0; k < num_samples; k++) {
        y_r[k] = H_lambda(r_raw[k], LAMBDA);
        y_i[k] = H_lambda(i_raw[k], LAMBDA);
    }

    // 3. Reconstruct
    // [중요] tid를 넘겨주어 core 함수도 쓰레드 전용 메모리를 쓰게 합니다.
    reconstruct_fast_core(y_r, rec_r, num_samples, LAMBDA, RECON_N, max_amp, tid);
    reconstruct_fast_core(y_i, rec_i, num_samples, LAMBDA, RECON_N, max_amp, tid);
    
    // ==================================================================
    // [핵심 수정] 2*Lambda Offset Correction (Anchor Point Fix)
    // CP/심볼 경계에서 발생하는 오프셋을 매번 0으로 강제 초기화
    // ==================================================================
    
    double two_lambda = 2.0 * LAMBDA;

    // (A) 실수부(Real) 보정
    // 복원된 첫 샘플과 입력 첫 샘플의 차이(Diff) 계산
    double diff_r = rec_r[0] - r_raw[0];
    
    // 차이가 2*Lambda의 정수배라면 그만큼 전체 이동
    // 예: diff가 0.1(2*L)이면 -> correction 0.1 -> 전체를 0.1 내림
    double correction_r = round(diff_r / two_lambda) * two_lambda;

    if (fabs(correction_r) > 1e-9) { 
        for (int k = 0; k < num_samples; k++) {
            rec_r[k] -= correction_r;
        }
    }

    // (B) 허수부(Imag) 보정
    double diff_i = rec_i[0] - i_raw[0];
    double correction_i = round(diff_i / two_lambda) * two_lambda;

    if (fabs(correction_i) > 1e-9) {
        for (int k = 0; k < num_samples; k++) {
            rec_i[k] -= correction_i;
        }
    }
    // ==================================================================
    
    // 4. Write Back (Saturation)
    for (int k = 0; k < num_samples; k++) {
        double clean_r = rec_r[k]; 
        double clean_i = rec_i[k];
        
        int32_t ir = (int32_t)(clean_r * 32768.0);
        int32_t ii = (int32_t)(clean_i * 32768.0);
        
        buffer[k].r = (int16_t)(ir > 32767 ? 32767 : (ir < -32768 ? -32768 : ir));
        buffer[k].i = (int16_t)(ii > 32767 ? 32767 : (ii < -32768 ? -32768 : ii));
    }
}

// ==================================================================================
// [New Function] Spatial Filtering (Interference Nulling)
// 입력: ru (OAI RU 구조체), start_idx (시작 인덱스), length (길이)
// 동작: ANT0과 ANT1의 데이터를 읽어서 간섭을 제거하고, 결과를 ANT0 버퍼에 덮어씁니다.
// ==================================================================================

void apply_spatial_filtering(RU_t* ru, int start_idx, int length) {
    
    // 1. 데이터 포인터 가져오기
    // OAI 메모리 구조에 맞춰 int32_t* -> c16_t* 캐스팅
    int32_t* ant0_ptr_32 = (int32_t*)ru->common.rxdata[0];
    int32_t* ant1_ptr_32 = (int32_t*)ru->common.rxdata[1];
    
    c16_t* rx0 = (c16_t*)&ant0_ptr_32[start_idx];
    c16_t* rx1 = (c16_t*)&ant1_ptr_32[start_idx];

    // 공분산 행렬 (Covariance Matrix) 요소
    double R00 = 0.0;   // |y0|^2
    // R11은 Nulling Weight 계산에 직접 안 쓰이므로 제거해도 되지만, 
    // 추후 디버깅(에너지 비교)을 위해 남겨두거나 삭제 가능. (여기선 삭제하여 최적화)
    // double R11 = 0.0; 
    
    double R01_r = 0.0; // Real(y0 * y1*)
    double R01_i = 0.0; // Imag(y0 * y1*)

    // -------------------------------------------------------
    // [Step 1] 공분산 행렬 계산 (Correlation)
    // -------------------------------------------------------
    for (int k = 0; k < length; k++) {
        double y0_r = (double)rx0[k].r;
        double y0_i = (double)rx0[k].i;
        double y1_r = (double)rx1[k].r;
        double y1_i = (double)rx1[k].i;

        // R00 = Sum(|y0|^2)
        R00 += (y0_r * y0_r + y0_i * y0_i);
        
        // R01 = Sum(y0 * conj(y1)) -> (y0_r + j y0_i)(y1_r - j y1_i)
        // Real: y0_r*y1_r + y0_i*y1_i
        // Imag: y0_i*y1_r - y0_r*y1_i
        R01_r += (y0_r * y1_r + y0_i * y1_i);
        R01_i += (y0_i * y1_r - y0_r * y1_i);
    }

    // -------------------------------------------------------
    // [Step 2] Nulling Weight 계산
    // 직교 벡터 w = [ -conj(R01), R00 ]
    // -------------------------------------------------------
    
    // w0 = -conj(R01) = -(R01_r - j R01_i) = -R01_r + j R01_i
    double w0_r = -R01_r;   
    double w0_i =  R01_i;   
    
    // w1 = R00 (실수)
    double w1_r = R00;         
    double w1_i = 0.0;

    // Weight 정규화 (Normalization)
    double norm_sq = w0_r*w0_r + w0_i*w0_i + w1_r*w1_r; // w1_i is 0
    double norm = sqrt(norm_sq);
    
    if (norm < 1e-9) norm = 1.0; // 0 나누기 방지

    double inv_norm = 1.0 / norm;
    w0_r *= inv_norm; 
    w0_i *= inv_norm;
    w1_r *= inv_norm; 
    // w1_i는 어차피 0이므로 연산 생략

    // -------------------------------------------------------
    // [Step 3] 필터링 적용 (Projection)
    // -------------------------------------------------------
    for (int k = 0; k < length; k++) {
        double y0_r = (double)rx0[k].r;
        double y0_i = (double)rx0[k].i;
        double y1_r = (double)rx1[k].r;
        double y1_i = (double)rx1[k].i;

        // Complex Multiply & Add: out = (w0 * y0) + (w1 * y1)
        // w1은 실수(Real)이므로 연산 간소화 가능
        
        // Term 1: w0 * y0
        double t1_r = w0_r * y0_r - w0_i * y0_i;
        double t1_i = w0_r * y0_i + w0_i * y0_r;

        // Term 2: w1 * y1 (w1_i == 0)
        double t2_r = w1_r * y1_r; 
        double t2_i = w1_r * y1_i;

        double out_r = t1_r + t2_r;
        double out_i = t1_i + t2_i;

        // Saturation Logic (int16 범위 클리핑)
        if (out_r > 32767.0) out_r = 32767.0;
        else if (out_r < -32768.0) out_r = -32768.0;

        if (out_i > 32767.0) out_i = 32767.0;
        else if (out_i < -32768.0) out_i = -32768.0;

        // 결과 저장 (ANT0 덮어쓰기)
        rx0[k].r = (int16_t)out_r;
        rx0[k].i = (int16_t)out_i;
        
        // ANT1 비우기 (Zeroing)
        rx1[k].r = 0;
        rx1[k].i = 0;
    }
}

// ==================================================================================
// [Phase Rotation / De-rotation Functions for Interference Shifting]
// ==================================================================================

// ==================================================================================
// [초고속 NCO (Look-Up Table) 기반 주파수 이동 모듈]
// ==================================================================================

#define NCO_LUT_SIZE 4096
float nco_cos_lut[NCO_LUT_SIZE];
float nco_sin_lut[NCO_LUT_SIZE];
int nco_initialized = 0;

// 시스템 시작 시 딱 한 번 호출하여 정답지(LUT)를 만듭니다.
void init_nco_lut(void) {
    if (nco_initialized) return;
    for (int i = 0; i < NCO_LUT_SIZE; i++) {
        double theta = 2.0 * M_PI * i / NCO_LUT_SIZE;
        nco_cos_lut[i] = (float)cos(theta);
        nco_sin_lut[i] = (float)sin(theta);
    }
    nco_initialized = 1;
}

// ----------------------------------------------------------------------------------
// 1. apply_rotation: 초고속 정방향 회전 (Spatial Filtering 이전)
// ----------------------------------------------------------------------------------
void apply_rotation(RU_t* ru, int start_idx, int length, uint64_t timestamp_rx, int sym_offset, double f_offset) {
    if (!nco_initialized) init_nco_lut();

    // 포인터에 restrict를 붙여 컴파일러의 AVX2 병렬화를 유도합니다.
    int32_t* restrict ant0_ptr_32 = (int32_t*)ru->common.rxdata[0];
    int32_t* restrict ant1_ptr_32 = (int32_t*)ru->common.rxdata[1];
    c16_t* restrict rx0 = (c16_t*)&ant0_ptr_32[start_idx];
    c16_t* restrict rx1 = (c16_t*)&ant1_ptr_32[start_idx];

    double sample_rate = ru->openair0_cfg.sample_rate; 
    
    // [핵심] 32-bit Phase Accumulator Step 계산
    // phase_step = (f_offset / fs) * 2^32
    // f_offset이 음수일 경우 2의 보수로 자동 변환되어 역회전합니다.
    uint32_t phase_step = (uint32_t)((f_offset / sample_rate) * 4294967296.0);

    // 절대 샘플 인덱스에 따른 현재 위상값 초기화
    uint64_t absolute_start_idx = timestamp_rx + sym_offset;
    uint32_t phase_acc = (uint32_t)(absolute_start_idx * phase_step);

    for (int k = 0; k < length; k++) {
        // 상위 12비트만 추출하여 LUT 인덱스(0~4095)로 사용
        uint32_t lut_idx = phase_acc >> 20;
        
        // 메모리에서 미리 계산된 cos, sin 값을 꺼내옴 (초고속)
        float curr_cos = nco_cos_lut[lut_idx];
        float curr_sin = nco_sin_lut[lut_idx];
        
        // --- Antenna 0 ---
        float i_val0 = (float)rx0[k].r;
        float q_val0 = (float)rx0[k].i;
        
        float out0_r = i_val0 * curr_cos - q_val0 * curr_sin;
        float out0_i = i_val0 * curr_sin + q_val0 * curr_cos;
        
        // --- Antenna 1 ---
        float i_val1 = (float)rx1[k].r;
        float q_val1 = (float)rx1[k].i;
        
        float out1_r = i_val1 * curr_cos - q_val1 * curr_sin;
        float out1_i = i_val1 * curr_sin + q_val1 * curr_cos;
        
        // Saturation (최적화를 위해 분기문 대신 3항 연산자 사용)
        out0_r = out0_r > 32767.0f ? 32767.0f : (out0_r < -32768.0f ? -32768.0f : out0_r);
        out0_i = out0_i > 32767.0f ? 32767.0f : (out0_i < -32768.0f ? -32768.0f : out0_i);
        out1_r = out1_r > 32767.0f ? 32767.0f : (out1_r < -32768.0f ? -32768.0f : out1_r);
        out1_i = out1_i > 32767.0f ? 32767.0f : (out1_i < -32768.0f ? -32768.0f : out1_i);

        rx0[k].r = (int16_t)out0_r; rx0[k].i = (int16_t)out0_i;
        rx1[k].r = (int16_t)out1_r; rx1[k].i = (int16_t)out1_i;

        // [핵심] 위상 업데이트 (무거운 곱셈 없이 오직 정수 덧셈 1번!)
        phase_acc += phase_step;
    }
}

// ----------------------------------------------------------------------------------
// 2. apply_derotation: 초고속 역방향 회전 (Spatial Filtering 이후)
// ----------------------------------------------------------------------------------
void apply_derotation(RU_t* ru, int start_idx, int length, uint64_t timestamp_rx, int sym_offset, double f_offset) {
    if (!nco_initialized) init_nco_lut();

    int32_t* restrict ant0_ptr_32 = (int32_t*)ru->common.rxdata[0];
    c16_t* restrict rx0 = (c16_t*)&ant0_ptr_32[start_idx];

    double sample_rate = ru->openair0_cfg.sample_rate; 
    
    // 역방향이므로 -f_offset을 줍니다.
    uint32_t phase_step = (uint32_t)((-f_offset / sample_rate) * 4294967296.0);

    uint64_t absolute_start_idx = timestamp_rx + sym_offset;
    uint32_t phase_acc = (uint32_t)(absolute_start_idx * phase_step);

    for (int k = 0; k < length; k++) {
        uint32_t lut_idx = phase_acc >> 20;
        
        float curr_cos = nco_cos_lut[lut_idx];
        float curr_sin = nco_sin_lut[lut_idx];
        
        float i_val0 = (float)rx0[k].r;
        float q_val0 = (float)rx0[k].i;
        
        float out0_r = i_val0 * curr_cos - q_val0 * curr_sin;
        float out0_i = i_val0 * curr_sin + q_val0 * curr_cos;
        
        out0_r = out0_r > 32767.0f ? 32767.0f : (out0_r < -32768.0f ? -32768.0f : out0_r);
        out0_i = out0_i > 32767.0f ? 32767.0f : (out0_i < -32768.0f ? -32768.0f : out0_i);

        rx0[k].r = (int16_t)out0_r;
        rx0[k].i = (int16_t)out0_i;

        phase_acc += phase_step;
    }
}
