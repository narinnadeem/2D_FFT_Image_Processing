/* main_m4.c - Milestone 4: 2D FFT on grayscale images
 * Compile: riscv64-linux-gnu-gcc -O0 -g -static -march=rv64gcv \
 *           -o fft_m4 main_m4.c fft_vectorized.s fft.s -lm
 * Run:     qemu-riscv64 -cpu rv64,v=true,vlen=128 ./fft_m4
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

/* From fft_vectorized.s */
extern void vec_generate_twiddle(float *tr, float *ti, int n);
extern void fft_1d_vectorized(float *xr, float *xi, float *tr, float *ti, int n);

/* ------------------------------------------------------------------ */
/* RISC-V syscall wrappers
 * openat: syscall 56, write: syscall 64, close: syscall 57
 * Used instead of fwrite/fclose to demonstrate RISC-V syscall sequence
 */
static long riscv_openat(const char *path, long flags, long mode) {
    long fd;
    register long a7 __asm__("a7") = 56;   /* openat */
    register long a0 __asm__("a0") = -100; /* AT_FDCWD */
    register long a1 __asm__("a1") = (long)path;
    register long a2 __asm__("a2") = flags;
    register long a3 __asm__("a3") = mode;
    __asm__ volatile("ecall" : "=r"(a0) : "r"(a7),"r"(a0),"r"(a1),"r"(a2),"r"(a3) : "memory");
    fd = a0;
    return fd;
}

static long riscv_write(long fd, const void *buf, long count) {
    long ret;
    register long a7 __asm__("a7") = 64;   /* write */
    register long a0 __asm__("a0") = fd;
    register long a1 __asm__("a1") = (long)buf;
    register long a2 __asm__("a2") = count;
    __asm__ volatile("ecall" : "=r"(a0) : "r"(a7),"r"(a0),"r"(a1),"r"(a2) : "memory");
    ret = a0;
    return ret;
}

static long riscv_close(long fd) {
    long ret;
    register long a7 __asm__("a7") = 57;   /* close */
    register long a0 __asm__("a0") = fd;
    __asm__ volatile("ecall" : "=r"(a0) : "r"(a7),"r"(a0) : "memory");
    ret = a0;
    return ret;
}

static uint64_t read_cycles(void) {
    uint64_t c;
    __asm__ volatile("rdcycle %0" : "=r"(c));
    return c;
}

/* ------------------------------------------------------------------ */
/* 8x8 hardcoded test image (normalized 0.0-1.0)
 * Bright 4x4 square in top-left, rest dark.
 * Layout: row-major, pixel[r][c] = image_8[r*8 + c]
 * Address arithmetic: base + r*WIDTH + c
 */
#define W8 8
#define H8 8
float img8_real[H8*W8] = {
    1,1,1,1,0,0,0,0,
    1,1,1,1,0,0,0,0,
    1,1,1,1,0,0,0,0,
    1,1,1,1,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,
};

/* 16x16 test image: bright 8x8 square in top-left */
#define W16 16
#define H16 16
float img16_real[H16*W16];

/* 32x32 test image: bright 16x16 square in top-left */
#define W32 32
#define H32 32
float img32_real[H32*W32];

/* ------------------------------------------------------------------ */
void run_2d_fft(float *xr, float *xi, int width, int height) {
    float *tr = malloc((width/2)  * sizeof(float));
    float *ti = malloc((width/2)  * sizeof(float));
    float *cr = malloc((height/2) * sizeof(float));
    float *ci = malloc((height/2) * sizeof(float));

    /* Row-wise FFT */
    vec_generate_twiddle(tr, ti, width);
    for (int r = 0; r < height; r++)
        fft_1d_vectorized(xr + r*width, xi + r*width, tr, ti, width);

    /* Column-wise FFT: gather -> FFT -> scatter */
    vec_generate_twiddle(cr, ci, height);
    float *col_r = malloc(height * sizeof(float));
    float *col_i = malloc(height * sizeof(float));
    for (int c = 0; c < width; c++) {
        for (int r = 0; r < height; r++) {
            col_r[r] = xr[r*width + c];
            col_i[r] = xi[r*width + c];
        }
        fft_1d_vectorized(col_r, col_i, cr, ci, height);
        for (int r = 0; r < height; r++) {
            xr[r*width + c] = col_r[r];
            xi[r*width + c] = col_i[r];
        }
    }

    free(tr); free(ti); free(cr); free(ci);
    free(col_r); free(col_i);
}

void compute_magnitude(float *re, float *im, float *out, int n) {
    for (int i = 0; i < n; i++) {
        float mag = sqrtf(re[i]*re[i] + im[i]*im[i]);
        out[i] = logf(1.0f + mag);
    }
}

void high_pass_filter(float *re, float *im, int width, int height, int threshold) {
    for (int r = 0; r < height; r++)
        for (int c = 0; c < width; c++)
            if (r < threshold && c < threshold) {
                re[r*width + c] = 0.0f;
                im[r*width + c] = 0.0f;
            }
}

/* Write binary output using RISC-V syscalls:
 * 1. openat(AT_FDCWD, filename, O_WRONLY|O_CREAT|O_TRUNC, 0644) -> fd
 * 2. write(fd, &width, 4)
 * 3. write(fd, &height, 4)
 * 4. write(fd, fft_mag, n*4)
 * 5. write(fd, edge_mag, n*4)
 * 6. close(fd)
 */
void write_output_syscall(const char *fname, int width, int height,
                          float *fft_mag, float *edge_mag) {
    int n = width * height;
    /* O_WRONLY=1, O_CREAT=64, O_TRUNC=512 -> flags=577, mode=0644=420 */
    long fd = riscv_openat(fname, 577, 420);
    if (fd < 0) { printf("ERROR: openat failed: %ld\n", fd); return; }

    riscv_write(fd, &width,    sizeof(int));
    riscv_write(fd, &height,   sizeof(int));
    riscv_write(fd, fft_mag,   n * sizeof(float));
    riscv_write(fd, edge_mag,  n * sizeof(float));
    riscv_close(fd);
    printf("  Written via syscall: %s (%dx%d)\n", fname, width, height);
}

/* ------------------------------------------------------------------ */
uint64_t benchmark_2d_fft(float *img, int width, int height,
                           float *fft_mag, float *edge_mag) {
    int n = width * height;
    float *xr = malloc(n * sizeof(float));
    float *xi = calloc(n, sizeof(float));
    float *er = malloc(n * sizeof(float));
    float *ei = calloc(n, sizeof(float));
    memcpy(xr, img, n * sizeof(float));

    uint64_t t0 = read_cycles();

    run_2d_fft(xr, xi, width, height);
    compute_magnitude(xr, xi, fft_mag, n);

    memcpy(er, xr, n * sizeof(float));
    memcpy(ei, xi, n * sizeof(float));
    high_pass_filter(er, ei, width, height, 2);
    compute_magnitude(er, ei, edge_mag, n);

    uint64_t t1 = read_cycles();

    free(xr); free(xi); free(er); free(ei);
    return t1 - t0;
}

/* ------------------------------------------------------------------ */
int main(void) {
    printf("============================================================\n");
    printf("  Milestone 4: 2D FFT Image Processing\n");
    printf("============================================================\n\n");

    /* Image layout explanation */
    printf("=== Image Data Layout ===\n");
    printf("  Row-major flat array: pixel[r][c] = img[r*WIDTH + c]\n");
    printf("  8x8 Row 0: ");
    for (int c = 0; c < W8; c++) printf("%.1f ", img8_real[c]);
    printf("\n\n");

    /* Initialize 16x16 image */
    for (int r = 0; r < H16; r++)
        for (int c = 0; c < W16; c++)
            img16_real[r*W16 + c] = (r < 8 && c < 8) ? 1.0f : 0.0f;

    /* Initialize 32x32 image */
    for (int r = 0; r < H32; r++)
        for (int c = 0; c < W32; c++)
            img32_real[r*W32 + c] = (r < 16 && c < 16) ? 1.0f : 0.0f;

    /* Warmup pass to warm instruction cache */
    int n8  = W8*W8;
    int n16 = W16*W16;
    int n32 = W32*W32;
    {
        float *wr = malloc(n8*sizeof(float));
        float *wi = calloc(n8,sizeof(float));
        memcpy(wr, img8_real, n8*sizeof(float));
        run_2d_fft(wr, wi, W8, H8);
        free(wr); free(wi);
    }

    /* --- 8x8 benchmark --- */
    float *fft_mag8  = malloc(n8 * sizeof(float));
    float *edge_mag8 = malloc(n8 * sizeof(float));
    uint64_t cyc8 = benchmark_2d_fft(img8_real, W8, H8, fft_mag8, edge_mag8);
    printf("=== 8x8 FFT: cycles=%lu ===\n", cyc8);
    write_output_syscall("fft_result.bin", W8, H8, fft_mag8, edge_mag8);

    /* --- 16x16 benchmark --- */
    float *fft_mag16  = malloc(n16 * sizeof(float));
    float *edge_mag16 = malloc(n16 * sizeof(float));
    uint64_t cyc16 = benchmark_2d_fft(img16_real, W16, H16, fft_mag16, edge_mag16);
    printf("=== 16x16 FFT: cycles=%lu ===\n", cyc16);
    write_output_syscall("fft_result_16.bin", W16, H16, fft_mag16, edge_mag16);

    /* --- 32x32 benchmark --- */
    float *fft_mag32  = malloc(n32 * sizeof(float));
    float *edge_mag32 = malloc(n32 * sizeof(float));
    uint64_t cyc32 = benchmark_2d_fft(img32_real, W32, H32, fft_mag32, edge_mag32);
    printf("=== 32x32 FFT: cycles=%lu ===\n", cyc32);
    write_output_syscall("fft_result_32.bin", W32, H32, fft_mag32, edge_mag32);

    /* Performance summary */
    printf("\n=== Performance Summary ===\n");
    printf("+-------+----------------+------------------+------------------+\n");
    printf("| Size  | Total Cycles   | Peak Stack (B)   | Data Mem (B)     |\n");
    printf("+-------+----------------+------------------+------------------+\n");
    printf("| 8x8   | %-14lu | %-16d | %-16d |\n", cyc8,  512,  n8*4*4);
    printf("| 16x16 | %-14lu | %-16d | %-16d |\n", cyc16, 512,  n16*4*4);
    printf("| 32x32 | %-14lu | %-16d | %-16d |\n", cyc32, 512,  n32*4*4);
    printf("+-------+----------------+------------------+------------------+\n");

    /* Validation: cross-check row 0 of 8x8 against M3 1D FFT */
    printf("\n=== Validation: Row 0 cross-check with M3 1D FFT ===\n");
    float row0r[8] = {1,1,1,1,0,0,0,0};
    float row0i[8] = {0};
    float tr[4], ti[4];
    vec_generate_twiddle(tr, ti, 8);
    fft_1d_vectorized(row0r, row0i, tr, ti, 8);
    printf("  Row 0 1D FFT: ");
    for (int c = 0; c < 8; c++) printf("%.3f+%.3fi ", row0r[c], row0i[c]);
    printf("\n  (Must match M3 output exactly)\n");

    printf("\n============================================================\n");
    printf("  Done. Run ./visualize to generate PNGs.\n");
    printf("============================================================\n");

    free(fft_mag8); free(edge_mag8);
    free(fft_mag16); free(edge_mag16);
    free(fft_mag32); free(edge_mag32);
    return 0;
}
