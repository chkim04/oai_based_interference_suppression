#ifndef RECONSTRUCTION_TEST_H
#define RECONSTRUCTION_TEST_H

#include "PHY/types.h"
#include "PHY/TOOLS/tools_defs.h"
#include <stdint.h>

/**
 * @brief Runs the Real-World Interference Reconstruction experiment.
 *
 * This function:
 * 1. Reads the real input buffer (PUSCH + Real Chirp).
 * 2. Measures the signal amplitude (beta_g).
 * 3. Emulates the Folding ADC (Modulo Operation) on the real data.
 * 4. Runs Algorithm 1 to unfold the signal.
 * 5. Removes the DC component (Blind Interference Suppression). ????????
 * 6. Overwrites the input buffer with the recovered PUSCH signal.
 *
 * @param buffer Pointer to the OAI c16_t buffer (Time Domain).
 * @param num_samples Number of samples to process.
 */
void run_real_interference_experiment(c16_t* buffer, int num_samples, uint32_t timestamp);

#endif // RECONSTRUCTION_TEST_H
