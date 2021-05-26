/**
 * @file fft.h
 *
 * @author Johannes Ehala, ProLab.
 * @license MIT
 *
 * Copyright ProLab, TTÜ. 26 Mai 2021
 */

#ifndef FFT_H_
#define FFT_H_

// Public functions
uint32_t calc_power (uint32_t num_samples);
void fft_complex (float *signal, uint32_t num_samples);

// Private functions
static void bit_alignment (float *data_startpoint, float *y_startpoint, uint32_t data_points);

#endif // FFT_H_ */
