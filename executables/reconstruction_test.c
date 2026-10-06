#include "reconstruction_test.h"
#include "common/utils/LOG/log.h"
#include "PHY/defs_common.h" 
#include <math.h>
#include <stdlib.h>
#include <string.h>

// --- CONFIGURATION ---
#define LAMBDA          0.0520       // ADC Threshold (Artificial Software Threshold)
#define RECON_N         1          // Difference Order (N=2 for Chirp)

// Scaling for q15 format
static const double q15_scale = 1.0 / 32768.0;

// --- HELPER FUNCTIONS ---

static double H_lambda(double input, double lambda) {
    double two_lambda = 2.0 * lambda;
    double frac_part = (input / two_lambda + 0.5) - floor(input / two_lambda + 0.5);
    return two_lambda * (frac_part - 0.5);
}

static double RD(double x, double lambda) {
    double two_lambda = 2.0 * lambda;
    return round(x / two_lambda) * two_lambda;
}

static void compute_diff(double* in, double* out, int size) {
    for (int i = 0; i < size - 1; i++) {
        out[i] = in[i+1] - in[i];
    }
    out[size - 1] = 0; 
}

static void cumsum(double* in, double* out, int size) {
    out[0] = 0; 
    for (int i = 0; i < size; i++) {
        out[i+1] = out[i] + in[i];
    }
}

// --- CORE ALGORITHM 1 (Same as before) ---
static void reconstruct_signal_core(double* y, double* recovered, int size, double lambda, int N, double beta_g) {
    int buf_size = size + 5; 
    double* ds  = (double*)calloc(buf_size, sizeof(double));
    double* rec = (double*)calloc(buf_size, sizeof(double));
    double* res = (double*)calloc(buf_size, sizeof(double));
    double* bl  = (double*)calloc(buf_size, sizeof(double));
    double* egN = (double*)calloc(buf_size, sizeof(double));
    double* s1  = (double*)calloc(buf_size, sizeof(double));
    double* s2  = (double*)calloc(buf_size, sizeof(double));
    double* tmp = (double*)calloc(buf_size, sizeof(double));

    // 1. Compute N-th Difference
    memcpy(ds, y, size * sizeof(double));
    int curr_len = size;
    for(int k=0; k<N; k++) {
        compute_diff(ds, tmp, curr_len);
        memcpy(ds, tmp, curr_len * sizeof(double));
        curr_len--; 
    }

    // 2. Modulo of Difference & Residual
    for(int i=0; i<curr_len; i++) {
        rec[i] = H_lambda(ds[i], lambda);
        res[i] = RD(rec[i] - ds[i], lambda);
    }

    // 3. Backward Reconstruction Loop
    memcpy(bl, res, curr_len * sizeof(double));
    memcpy(egN, res, curr_len * sizeof(double));
    int bl_len = curr_len;
    int J = (int)ceil(6.0 * beta_g / lambda);

    for (int k = 1; k < N; k++) {
        cumsum(egN, s1, bl_len);      
        cumsum(s1, s2, bl_len + 1);   
        
        int idx_1 = 1; 
        int idx_J = J + 1;
        if (idx_J >= bl_len + 2) idx_J = bl_len + 1; 

        double term1 = -s2[idx_J]; 
        double term2 = s2[idx_1]; 
        double Kn = floor((term1 + term2) / (12.0 * beta_g) + 0.5);

        cumsum(bl, tmp, bl_len);
        bl_len++; 
        
        double correction = 2.0 * Kn * lambda;
        for(int i=0; i<bl_len; i++) {
            bl[i] = RD(tmp[i], lambda) + correction;
        }
        memcpy(egN, bl, bl_len * sizeof(double));
    }

    // 4. Final Integration
    cumsum(bl, tmp, bl_len);
    for(int i=0; i<size; i++) {
        recovered[i] = tmp[i] + y[i];
    }

    free(ds); free(rec); free(res); free(bl); free(egN); free(s1); free(s2); free(tmp);
}

// --- MAIN ENTRY POINT (REAL WORLD INTERFERENCE) ---
void run_real_interference_experiment(c16_t* buffer, int num_samples, uint32_t timestamp) {
    
    // Allocate Buffers
    double* real_raw = (double*)malloc(num_samples * sizeof(double));
    double* imag_raw = (double*)malloc(num_samples * sizeof(double));
    double* y_r      = (double*)malloc(num_samples * sizeof(double));
    double* y_i      = (double*)malloc(num_samples * sizeof(double));
    double* rec_r    = (double*)malloc(num_samples * sizeof(double));
    double* rec_i    = (double*)malloc(num_samples * sizeof(double));

    // 1. Load Real Hardware Signal & Measure Max Amplitude
    double max_amp_measured = 0.0;

    for (int i = 0; i < num_samples; i++) {
        // Convert q15 to double
        real_raw[i] = (double)buffer[i].r * q15_scale;
        imag_raw[i] = (double)buffer[i].i * q15_scale;

        // Track Peak Amplitude (Beta_g)
        // This replaces the manual 'CHIRP_AMP' setting
        if(fabs(real_raw[i]) > max_amp_measured) max_amp_measured = fabs(real_raw[i]);
        if(fabs(imag_raw[i]) > max_amp_measured) max_amp_measured = fabs(imag_raw[i]);
    }

    // SAFETY CHECK: If signal is too weak, algorithm might become unstable
    // Ensure beta_g is at least slightly larger than lambda
    if (max_amp_measured < LAMBDA) max_amp_measured = LAMBDA * 1.1;

    // 2. Emulate Folding (The "Virtual ADC")
    // We artificially fold the strong real signal into [-lambda, lambda]
    for (int i = 0; i < num_samples; i++) {
        y_r[i] = H_lambda(real_raw[i], LAMBDA);
        y_i[i] = H_lambda(imag_raw[i], LAMBDA);
        
        // --- DEBUG: WRITE BACK FOLDED DATA TO OAI BUFFER ---
    	// Do this temporarily to visualize the "Folded" state in T-view
    	// This overwrites the OAI buffer with the FOLDED signal.
    	//buffer[i].r = (int16_t)(y_r[i] * 32768.0);
    	//buffer[i].i = (int16_t)(y_i[i] * 32768.0);
    }
    
    /*
    // We create a temporary buffer to pack the data for the T-Tracer
    // without touching the main 'buffer' or 'y_r' arrays used for processing.
    int32_t* viz_buffer = (int32_t*)malloc(num_samples * sizeof(int32_t));
    
    if (viz_buffer) {
        for (int k = 0; k < num_samples; k++) {
            // Scale double back to int16 range for visualization
            int16_t i_val = (int16_t)(y_r[k] * 32768.0);
            int16_t q_val = (int16_t)(y_i[k] * 32768.0);
            
            // Pack into OAI's 32-bit IQ format
            viz_buffer[k] = (i_val & 0xFFFF) | (q_val << 16);
        }

        // SEND TO T-TRACER
        // ID: T_USRP_RX_ANT0 (As requested)
        // Format: int(module/ant), int(timestamp), buffer
        T(T_USRP_RX_ANT0, 
          T_INT(timestamp),      // Correct Time Axis
          T_BUFFER(viz_buffer, num_samples * 4)); // 4 bytes per sample

        free(viz_buffer);
    }
    */
    
    
    //3. Reconstruct (Algorithm 1)
    reconstruct_signal_core(y_r, rec_r, num_samples, LAMBDA, RECON_N, max_amp_measured);
    reconstruct_signal_core(y_i, rec_i, num_samples, LAMBDA, RECON_N, max_amp_measured);
    
    /*
    // 4. Blind Interference Suppression (DC Subtraction)
    // Since we don't know the exact chirp phase to subtract perfectly,
    // we subtract the Mean (DC) which contains the bulk of the low-freq chirp.
    double sum_r = 0, sum_i = 0;
    for(int i=0; i<num_samples; i++) { sum_r += rec_r[i]; sum_i += rec_i[i]; }
    double dc_r = sum_r / num_samples;
    double dc_i = sum_i / num_samples;

    // 5. Write Back to OAI Buffer
    for (int i = 0; i < num_samples; i++) {
        // Remove DC (The Real Chirp)
        double clean_r = rec_r[i] - dc_r;
        double clean_i = rec_i[i] - dc_i;
        
        // Scale back to int16
        // Note: We use saturation logic just in case reconstruction spiked
        int32_t tmp_r = (int32_t)(clean_r * 32768.0);
        int32_t tmp_i = (int32_t)(clean_i * 32768.0);
        
        buffer[i].r = (int16_t)(tmp_r > 32767 ? 32767 : (tmp_r < -32768 ? -32768 : tmp_r));
        buffer[i].i = (int16_t)(tmp_i > 32767 ? 32767 : (tmp_i < -32768 ? -32768 : tmp_i));
    }
    */
    
    free(real_raw); free(imag_raw); 
    free(y_r); free(y_i); 
    free(rec_r); free(rec_i);
}
