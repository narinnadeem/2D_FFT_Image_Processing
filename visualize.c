/* visualize.c - Host-side PNG generator for Milestone 4
 * Reads fft_result.bin (and _16, _32 variants) and produces output PNGs.
 * Compile: gcc visualize.c -o visualize -lm
 * Run:     ./visualize
 */

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

static void normalize_to_gray(float *data, unsigned char *out, int n) {
    float mn = data[0], mx = data[0];
    for (int i = 1; i < n; i++) {
        if (data[i] < mn) mn = data[i];
        if (data[i] > mx) mx = data[i];
    }
    float range = mx - mn;
    if (range < 1e-10f) range = 1.0f;
    for (int i = 0; i < n; i++) {
        float norm = (data[i] - mn) / range;
        out[i] = (unsigned char)(norm * 255.0f + 0.5f);
    }
}

static void fftshift(float *data, float *out, int width, int height) {
    int hw = width/2, hh = height/2;
    for (int r = 0; r < height; r++)
        for (int c = 0; c < width; c++) {
            int nr = (r + hh) % height;
            int nc = (c + hw) % width;
            out[nr*width + nc] = data[r*width + c];
        }
}

int main(void) {
    const char *bins[]  = {"fft_result.bin", "fft_result_16.bin", "fft_result_32.bin"};
    const char *ffts[]  = {"fft_output.png", "fft_output_16.png", "fft_output_32.png"};
    const char *edges[] = {"edges_output.png", "edges_output_16.png", "edges_output_32.png"};

    for (int b = 0; b < 3; b++) {
        FILE *f = fopen(bins[b], "rb");
        if (!f) { printf("Skipping %s (not found)\n", bins[b]); continue; }

        int width, height;
        fread(&width,  sizeof(int), 1, f);
        fread(&height, sizeof(int), 1, f);
        int n = width * height;
        printf("Reading %s: %dx%d\n", bins[b], width, height);

        float *fft_mag  = malloc(n * sizeof(float));
        float *edge_mag = malloc(n * sizeof(float));
        fread(fft_mag,  sizeof(float), n, f);
        fread(edge_mag, sizeof(float), n, f);
        fclose(f);

        float *fft_shifted = malloc(n * sizeof(float));
        fftshift(fft_mag, fft_shifted, width, height);

        unsigned char *fft_gray  = malloc(n);
        unsigned char *edge_gray = malloc(n);
        normalize_to_gray(fft_shifted, fft_gray,  n);
        normalize_to_gray(edge_mag,    edge_gray, n);

        stbi_write_png(ffts[b],  width, height, 1, fft_gray,  width);
        stbi_write_png(edges[b], width, height, 1, edge_gray, width);
        printf("  Written: %s and %s\n", ffts[b], edges[b]);

        free(fft_mag); free(edge_mag); free(fft_shifted);
        free(fft_gray); free(edge_gray);
    }
    return 0;
}
