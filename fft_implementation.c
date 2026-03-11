#define PI 3.14159265358979323846f
#include "fft_implementation.h"
#include "math.h"
#include <stdio.h>
#include <stdbool.h>  
#include <stdlib.h>

unsigned int reverse_bits(unsigned int x, int log_n) {
    unsigned int reversed = 0;
    for (int i = 0; i < log_n; i++) {
        // Shift reversed left and bring in the LSB of x
        reversed = (reversed << 1) | (x & 1);
        x >>= 1; // Shift x right to get the next bit
    }
    return reversed;
}

void bit_reverse_array(float* x_real, float* x_imag, int n) {
    int log_n = log2_int(n); 
    for (int i = 0; i < n; i++) {
        unsigned int j = reverse_bits(i, log_n);
        
        // Only swap if i < j to avoid swapping elements back to original positions
        if (i < j) {
            // Swap real parts
            float temp_real = x_real[i];
            x_real[i] = x_real[j];
            x_real[j] = temp_real;
            
            // Swap imaginary parts
            float temp_imag = x_imag[i];
            x_imag[i] = x_imag[j];
            x_imag[j] = temp_imag;
        }
    }
}

void generate_twiddle_factors(float* twiddle_real, float* twiddle_imag, int n) {
    // need n/2 twiddle factors for an FFT of size n 
    for (int k = 0; k < n / 2; k++) {
        float angle = -2.0f * PI * k / n;
        twiddle_real[k] = mycos(angle); 
        twiddle_imag[k] = mysin(angle); 
    }
}

void butterfly_iterative(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n) {
    int stages = log2_int(n); 

    for (int stage = 1; stage <= stages; stage++) {
        int m = 1 << stage;           
        int half_m = m / 2;           
        int twiddle_step = n / m;     

        for (int k = 0; k < n; k += m) {
            for (int j = 0; j < half_m; j++) {
                int top = k + j;
                int bottom = k + j + half_m;

                float wr = twiddle_real[j * twiddle_step];
                float wi = twiddle_imag[j * twiddle_step];

                float tr = wr * x_real[bottom] - wi * x_imag[bottom];
                float ti = wr * x_imag[bottom] + wi * x_real[bottom];

                x_real[bottom] = x_real[top] - tr;
                x_imag[bottom] = x_imag[top] - ti;
                x_real[top] = x_real[top] + tr;
                x_imag[top] = x_imag[top] + ti;
            }
        }
    }
}

void fft_1d_iterative(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n) {
    // 1. Reorder the array into bit-reversed order 
    bit_reverse_array(x_real, x_imag, n);
    // 2. Perform the actual FFT computation
    butterfly_iterative(x_real, x_imag, twiddle_real, twiddle_imag, n);
}

void butterfly_recursive(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n, int total_n) {
    if (n <= 1) return;
    int half = n / 2;

    butterfly_recursive(x_real, x_imag, twiddle_real, twiddle_imag, half, total_n);
    butterfly_recursive(x_real + half, x_imag + half, twiddle_real, twiddle_imag, half, total_n);

    for (int j = 0; j < half; j++) {
        int t_idx = j * (total_n / n);
        float wr = twiddle_real[t_idx];
        float wi = twiddle_imag[t_idx];

        float tr = (x_real[j + half] * wr) - (x_imag[j + half] * wi);
        float ti = (x_real[j + half] * wi) + (x_imag[j + half] * wr);

        x_real[j + half] = x_real[j] - tr;
        x_imag[j + half] = x_imag[j] - ti;
        x_real[j] = x_real[j] + tr;
        x_imag[j] = x_imag[j] + ti;
    }
}

void fft_1d_recursive(float* x_real, float* x_imag, float* twiddle_real, float* twiddle_imag, int n) {
    bit_reverse_array(x_real, x_imag, n); 
    butterfly_recursive(x_real, x_imag, twiddle_real, twiddle_imag, n, n); 
}

void transpose(float** matrix, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = i + 1; j < cols; j++) {
            float temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

void fft_2d(float** x_real, float** x_imag, float* twiddle_real, float* twiddle_imag, int rows, int cols, bool test, bool inverse) {
    // 1. Generate twiddles for the current size (N = cols)
    generate_twiddle_factors(twiddle_real, twiddle_imag, cols);

    // 2. Apply 1D FFT to each row
    for (int i = 0; i < rows; i++) {
        fft_1d_iterative(x_real[i], x_imag[i], twiddle_real, twiddle_imag, cols);
    }

    // 3. Transpose the matrix 
    transpose(x_real, rows, cols);
    transpose(x_imag, rows, cols);
    
    // 4. Apply 1D FFT to each column
    for (int i = 0; i < rows; i++) {
        fft_1d_iterative(x_real[i], x_imag[i], twiddle_real, twiddle_imag, cols);
    }

    // 5. Transpose back to original orientation 
    transpose(x_real, rows, cols);
    transpose(x_imag, rows, cols);
}

void ifft_2d(float** x_real, float** x_imag, float* twiddle_real, float* twiddle_imag, int rows, int cols, bool test) {
    // 1. Generate inverse twiddle factors (imaginary signs flipped) 
    for (int k = 0; k < cols / 2; k++) {
        float angle = 2.0f * PI * k / cols; // Positive angle for inverse
        twiddle_real[k] = mycos(angle);
        twiddle_imag[k] = mysin(angle);
    }

    // 2. Call the 2D FFT with the inverse twiddles
    fft_2d(x_real, x_imag, twiddle_real, twiddle_imag, rows, cols, test, true);

    // 3. Normalize: divide every element by (rows * cols) 
    float num_elements = (float)(rows * cols);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            x_real[i][j] /= num_elements;
            x_imag[i][j] /= num_elements;
        }
    }
}

void create_highpass_filter(float** filter, int rows, int cols, float cutoff) {
    float d0_sq = 2.0f * cutoff * cutoff;
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Handle wrap-around distance (centering the filter)
            float di = (i > rows / 2) ? (float)(i - rows) : (float)i;
            float dj = (j > cols / 2) ? (float)(j - cols) : (float)j;
            
            // Calculate distance D(i,j) = sqrt(di^2 + dj^2) 
            float d_sq = (di * di) + (dj * dj);
            
            // H(i,j) = 1 - exp(-D^2 / 2D0^2) 
            filter[i][j] = 1.0f - myexp(-d_sq / d0_sq);
        }
    }
}

void edge_detection(float** x_real, float** x_imag, float* twiddle_real, float* twiddle_imag, float** filter, int rows, int cols, float cutoff) {
    // 1. Create the highpass filter 
    create_highpass_filter(filter, rows, cols, cutoff);
    // 2. Multiply FFT result by the filter 
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            x_real[i][j] *= filter[i][j];
            x_imag[i][j] *= filter[i][j];
        }
    }
    // 3. Inverse FFT to convert back to spatial domain 
    ifft_2d(x_real, x_imag, twiddle_real, twiddle_imag, rows, cols, false);
}
