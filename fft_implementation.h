#include <stdio.h>
#include <stdbool.h>  
#include <stdlib.h>

unsigned int reverse_bits(unsigned int x, int log_n);
void bit_reverse_array(float* x_real, float* x_imag, int n);
void generate_twiddle_factors(float* twiddle_real, float* twiddle_imag, int n);
void butterfly_iterative(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n);
void fft_1d_iterative(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n);
void butterfly_recursive(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n, int total_n);void fft_1d_recursive(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n);
void transpose(float** matrix, int rows, int cols);
void fft_2d(float** x_real, float** x_imag, float* twiddle_real, float* twiddle_imag, int rows, int cols, bool test, bool inverse);
void ifft_2d(float** x_real, float** x_imag, float* twiddle_real, float* twiddle_imag, int rows, int cols, bool test);
void create_highpass_filter(float** filter, int rows, int cols, float cutoff);
void edge_detection(float** x_real, float** x_imag, float* twiddle_real, float* twiddle_imag, float** filter, int rows, int cols, float cutoff);
