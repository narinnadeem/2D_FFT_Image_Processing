#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h> 
#include "stb_image.h"
#include "stb_image_write.h"
#include "math.h"
#include "fft_implementation.h"

int main() {
    int width, height, channels;

    // 1. Load the image
    unsigned char *img = stbi_load("OIP.png", &width, &height, &channels, 1); 
    if (img == NULL) {
        printf("Failed to load image\n"); 
        return 1;
    }

    // 2. Memory Allocation for 2D Arrays
    float **real = (float **)malloc(height * sizeof(float *));
    float **imag = (float **)malloc(height * sizeof(float *));
    for (int i = 0; i < height; i++) {
        real[i] = (float *)malloc(width * sizeof(float));
        imag[i] = (float *)malloc(width * sizeof(float));
    }

    // 3. Convert uint8 [0-255] to float [0.0-1.0]
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            real[i][j] = img[i * width + j] / 255.0f;
            imag[i][j] = 0.0f;
        }
    }

    // 4. Initialize Twiddles
    float *twiddle_real = (float *)malloc((width / 2) * sizeof(float));
    float *twiddle_imag = (float *)malloc((width / 2) * sizeof(float));
    generate_twiddle_factors(twiddle_real, twiddle_imag, width);

    // 5. Perform 2D FFT
    fft_2d(real, imag, twiddle_real, twiddle_imag, height, width, false, false);

    // We save this now because edge_detection will change these arrays
    unsigned char *output_data = (unsigned char *)malloc(width * height);
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            float mag = sqrt(real[i][j]*real[i][j] + imag[i][j]*imag[i][j]);
            // Log scaling: helps visualize intensities
            float scaled = 100.0f * logf(1.0f + mag); 
            if (scaled > 255.0f) scaled = 255.0f;
            output_data[i * width + j] = (unsigned char)scaled;
        }
    }
    stbi_write_png("fft_output.png", width, height, 1, output_data, width);

    // --- EDGE DETECTION ---
    // 1. Allocate memory for the filter
    float **filter = (float **)malloc(height * sizeof(float *));
    for (int i = 0; i < height; i++) {
        filter[i] = (float *)malloc(width * sizeof(float));
    }

    // 2. Create the Gaussian High-Pass Filter
    create_highpass_filter(filter, height, width, 50.0f);

    // 3. Apply the filter and perform Inverse FFT
edge_detection(real, imag, twiddle_real, twiddle_imag, filter, height, width, 1.5f);

    // 4. Save the Resulting Edges
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            // After IFFT, real part contains the spatial image. 
            // We scale it up because edges on 8x8 can be very faint.
            float val = real[i][j] * 255.0f; 
            if (val > 255.0f) val = 255.0f;
            if (val < 0.0f) val = 0.0f;
            output_data[i * width + j] = (unsigned char)val;
        }
    }
    stbi_write_png("edges_output.png", width, height, 1, output_data, width);
    printf("FFT Processing Complete. Created fft_output.png and edges_output.png\n"); 

    // 7. Cleanup
    for (int i = 0; i < height; i++) {
        free(filter[i]);
    }
    free(filter);
    free(img);
    free(output_data);
    free(twiddle_real);
    free(twiddle_imag);
    for (int i = 0; i < height; i++) {
        free(real[i]);
        free(imag[i]);
    }
    free(real);
    free(imag);

    return 0;
}