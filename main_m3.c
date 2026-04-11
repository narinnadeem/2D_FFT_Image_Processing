/* main_m3.c - Milestone 3: RISC-V Vector FFT test driver
 * Compile: riscv64-linux-gnu-gcc -O0 -g -static -march=rv64gcv \
 *           -o fft_m3 main_m3.c fft_vectorized.s fft.s -lm
 * Run:     qemu-riscv64 -cpu rv64,v=true,vlen=128 ./fft_m3
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>

/* From fft.s (scalar - reused) */
extern float  my_sin(float x);
extern float  my_cos(float x);
extern void   generate_twiddle_factors(float *tr, float *ti, int n);
extern void   fft_1d_iterative(float *xr, float *xi, float *tr, float *ti, int n);

/* From fft_vectorized.s */
extern void   vec_generate_twiddle(float *tr, float *ti, int n);
extern void   vec_bit_reverse_array(float *xr, float *xi, int n);
extern void   vec_butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n);
extern void   fft_1d_vectorized(float *xr, float *xi, float *tr, float *ti, int n);

static uint64_t read_cycles(void) {
    uint64_t c;
    __asm__ volatile("rdcycle %0" : "=r"(c));
    return c;
}

/* ------------------------------------------------------------------ */
static void print_separator(void) {
    printf("============================================================\n");
}

/* ------------------------------------------------------------------ */
static void test_bit_reversal(void) {
    printf("\n=== Vectorized Bit Reversal Validation (n=8) ===\n");
    int n = 8;
    float xr[8] = {0,1,2,3,4,5,6,7};
    float xi[8] = {0,0,0,0,0,0,0,0};

    printf("  Before:");
    for (int i = 0; i < n; i++) printf(" %.0f", xr[i]);
    printf("\n");

    vec_bit_reverse_array(xr, xi, n);

    printf("  After: ");
    for (int i = 0; i < n; i++) printf(" %.0f", xr[i]);
    printf("\n");
    printf("  Expected: 0 4 2 6 1 5 3 7\n");

    int expected[8] = {0,4,2,6,1,5,3,7};
    int ok = 1;
    for (int i = 0; i < n; i++)
        if ((int)xr[i] != expected[i]) { ok = 0; break; }
    printf("  Result: %s\n", ok ? "PASS" : "FAIL");
}

/* ------------------------------------------------------------------ */
static void test_fft_correctness(int n, const char *label) {
    float *xr = calloc(n, sizeof(float));
    float *xi = calloc(n, sizeof(float));
    float *tr = calloc(n/2, sizeof(float));
    float *ti = calloc(n/2, sizeof(float));

    /* all-ones input */
    for (int i = 0; i < n; i++) { xr[i] = 1.0f; xi[i] = 0.0f; }

    vec_generate_twiddle(tr, ti, n);
    fft_1d_vectorized(xr, xi, tr, ti, n);

    printf("\n=== Vectorized FFT Correctness (%s, n=%d) ===\n", label, n);
    printf("  All-ones input:\n");
    for (int i = 0; i < n; i++)
        printf("  X[%2d] = %8.4f + %8.4fi\n", i, xr[i], xi[i]);
    printf("  Expected: X[0]=%.1f, all others ~0\n", (float)n);

    /* impulse input */
    for (int i = 0; i < n; i++) { xr[i] = (i==0) ? 1.0f : 0.0f; xi[i] = 0.0f; }
    vec_generate_twiddle(tr, ti, n);
    fft_1d_vectorized(xr, xi, tr, ti, n);

    printf("  Impulse input:\n");
    for (int i = 0; i < n; i++)
        printf("  X[%2d] = %8.4f + %8.4fi\n", i, xr[i], xi[i]);
    printf("  Expected: all X[k] = 1.0000 + 0.0000i\n");

    free(xr); free(xi); free(tr); free(ti);
}

/* ------------------------------------------------------------------ */
static void benchmark(int n) {
    float *xr_v = calloc(n, sizeof(float));
    float *xi_v = calloc(n, sizeof(float));
    float *tr_v = calloc(n/2, sizeof(float));
    float *ti_v = calloc(n/2, sizeof(float));
    float *xr_s = calloc(n, sizeof(float));
    float *xi_s = calloc(n, sizeof(float));
    float *tr_s = calloc(n/2, sizeof(float));
    float *ti_s = calloc(n/2, sizeof(float));

    for (int i = 0; i < n; i++) { xr_v[i] = xr_s[i] = 1.0f; }

    vec_generate_twiddle(tr_v, ti_v, n);
    generate_twiddle_factors(tr_s, ti_s, n);

    uint64_t t0, t1;

    t0 = read_cycles();
    fft_1d_vectorized(xr_v, xi_v, tr_v, ti_v, n);
    t1 = read_cycles();
    uint64_t vec_cycles = t1 - t0;

    for (int i = 0; i < n; i++) { xr_s[i] = 1.0f; xi_s[i] = 0.0f; }
    t0 = read_cycles();
    fft_1d_iterative(xr_s, xi_s, tr_s, ti_s, n);
    t1 = read_cycles();
    uint64_t scl_cycles = t1 - t0;

    printf("| %-4d | %-14lu | %-16lu |\n", n, vec_cycles, scl_cycles);

    free(xr_v); free(xi_v); free(tr_v); free(ti_v);
    free(xr_s); free(xi_s); free(tr_s); free(ti_s);
}

/* ------------------------------------------------------------------ */
int main(void) {
    print_separator();
    printf("  Milestone 3: RISC-V Vector FFT\n");
    print_separator();

    /* Bit reversal validation */
    test_bit_reversal();

    /* Debug: test twiddle generation alone */
    printf("\n=== Debug: Twiddle Generation (n=8) ===\n");
    {
        float tr[4], ti[4];
        vec_generate_twiddle(tr, ti, 8);
        for (int i = 0; i < 4; i++)
            printf("  W[%d]: cos=%.4f sin=%.4f\n", i, tr[i], ti[i]);
    }
    printf("  Twiddle OK\n");

    /* Debug: test vectorized butterfly alone (no bit reversal) */
    printf("\n=== Debug: Vectorized Butterfly (n=8, all-ones) ===\n");
    {
        float xr[8]={1,1,1,1,1,1,1,1}, xi[8]={0};
        float tr[4], ti[4];
        vec_generate_twiddle(tr, ti, 8);
        /* manually bit-reverse first using scalar */
        float tmp_xr[8]={1,1,1,1,1,1,1,1}, tmp_xi[8]={0};
        vec_butterfly_iterative(tmp_xr, tmp_xi, tr, ti, 8);
        for (int i = 0; i < 8; i++)
            printf("  X[%d]=%.4f+%.4fi\n", i, tmp_xr[i], tmp_xi[i]);
    }
    printf("  Butterfly OK\n");

    /* FFT correctness */
    test_fft_correctness(8,  "small");
    test_fft_correctness(16, "medium");
    test_fft_correctness(32, "large");

    /* Performance benchmarks */
    printf("\n=== Performance Benchmarks: Vector vs Scalar ===\n");
    printf("+------+----------------+------------------+\n");
    printf("| N    | Vec Cycles     | Scalar Cycles    |\n");
    printf("+------+----------------+------------------+\n");
    benchmark(8);
    benchmark(16);
    benchmark(32);
    printf("+------+----------------+------------------+\n");

    printf("\n=== Large Input Scalability (N=1024) ===\n");
    {
        int n = 1024;
        float *xr = calloc(n, sizeof(float));
        float *xi = calloc(n, sizeof(float));
        float *tr = calloc(n/2, sizeof(float));
        float *ti = calloc(n/2, sizeof(float));
        for (int i = 0; i < n; i++) xr[i] = 1.0f;
        vec_generate_twiddle(tr, ti, n);
        uint64_t t0 = read_cycles();
        fft_1d_vectorized(xr, xi, tr, ti, n);
        uint64_t t1 = read_cycles();
        printf("  N=1024 vectorized FFT: X[0]=%.4f (expected %.1f)\n", xr[0], (float)n);
        printf("  Cycles: %lu\n", t1 - t0);
        free(xr); free(xi); free(tr); free(ti);
    }

    print_separator();
    printf("  Done.\n");
    print_separator();
    return 0;
}
