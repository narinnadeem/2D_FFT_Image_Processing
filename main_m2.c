/*
 * main.c  –  Milestone 2 driver
 * Calls RISC-V assembly FFT functions, validates output, measures performance.
 * Compile with:
 *   riscv64-linux-gnu-gcc -O0 -o fft_m2 main.c fft.s -lm
 *   qemu-riscv64 ./fft_m2
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/* ---------- Declarations of assembly functions ---------- */
extern float  my_sin(float x);
extern float  my_cos(float x);
extern void   generate_twiddle_factors(float *tr, float *ti, int n);
extern unsigned int reverse_bits(unsigned int x, int log_n);
extern void   bit_reverse_array(float *xr, float *xi, int n);
extern void   butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n);
extern void   fft_1d_iterative(float *xr, float *xi, float *tr, float *ti, int n);
extern void   butterfly_recursive(float *xr, float *xi, float *tr, float *ti, int n, int total_n);
extern void   fft_1d_recursive(float *xr, float *xi, float *tr, float *ti, int n);
extern int    log2_int(int n);

/* ---------- Cycle counter (RISC-V rdcycle) ---------- */
static inline unsigned long read_cycles(void) {
    unsigned long c;
    __asm__ volatile ("rdcycle %0" : "=r"(c));
    return c;
}

/* ---------- Helpers ---------- */
#define MAX_N 32

static void copy_input(float *dst_r, float *dst_i, const float *src, int n) {
    for (int i = 0; i < n; i++) { dst_r[i] = src[i]; dst_i[i] = 0.0f; }
}

static void print_fft(const char *label, const float *xr, const float *xi, int n) {
    printf("\n%s (n=%d):\n", label, n);
    for (int i = 0; i < n; i++) {
        printf("  X[%2d] = %8.4f + %8.4fi\n", i, xr[i], xi[i]);
    }
}

/* ---------- Math validation ---------- */
static void test_math(void) {
    printf("=== Math Function Validation ===\n");
    float angles[] = {0.0f, 0.7854f, 1.5708f, 3.1416f, -1.5708f};
    const char *names[] = {"0", "PI/4", "PI/2", "PI", "-PI/2"};
    for (int i = 0; i < 5; i++) {
        float s_asm = my_sin(angles[i]);
        float c_asm = my_cos(angles[i]);
        float s_ref = sinf(angles[i]);
        float c_ref = cosf(angles[i]);
        printf("  angle=%-6s  sin: asm=%7.4f ref=%7.4f | cos: asm=%7.4f ref=%7.4f\n",
               names[i], s_asm, s_ref, c_asm, c_ref);
    }
}

/* ---------- Bit-reversal validation ---------- */
static void test_bit_reversal(void) {
    printf("\n=== Bit Reversal Validation (n=8, log_n=3) ===\n");
    printf("  Index | Binary | Bit-Rev Binary | Bit-Rev Index | Assembly\n");
    printf("  ------|--------|----------------|---------------|----------\n");
    int expected[] = {0, 4, 2, 6, 1, 5, 3, 7};
    for (int i = 0; i < 8; i++) {
        unsigned int r = reverse_bits((unsigned int)i, 3);
        printf("  %5d |  %03d   |      %03d       |      %d        | %u %s\n",
               i, i, expected[i], expected[i], r,
               (r == (unsigned)expected[i]) ? "OK" : "FAIL");
    }
}

/* ---------- Array-reorder validation ---------- */
static void test_array_reorder(void) {
    printf("\n=== Array Reorder Validation (n=8) ===\n");
    float xr[8] = {0,1,2,3,4,5,6,7};
    float xi[8] = {0};
    printf("  Before: ");
    for (int i = 0; i < 8; i++) printf("%.0f ", xr[i]);
    printf("\n");
    bit_reverse_array(xr, xi, 8);
    printf("  After:  ");
    for (int i = 0; i < 8; i++) printf("%.0f ", xr[i]);
    printf("\n  Expected: 0 4 2 6 1 5 3 7\n");
}

/* ---------- FFT validation ---------- */
static void test_fft(int n) {
    /* Simple DC signal: all ones → X[0]=n, rest=0 */
    float input[MAX_N];
    for (int i = 0; i < n; i++) input[i] = 1.0f;

    float tr[MAX_N/2], ti[MAX_N/2];
    generate_twiddle_factors(tr, ti, n);

    /* Iterative */
    float xr_it[MAX_N], xi_it[MAX_N];
    copy_input(xr_it, xi_it, input, n);
    fft_1d_iterative(xr_it, xi_it, tr, ti, n);
    print_fft("Iterative FFT (all-ones input)", xr_it, xi_it, n);

    /* Recursive – need fresh twiddles (they're not consumed) */
    generate_twiddle_factors(tr, ti, n);
    float xr_rc[MAX_N], xi_rc[MAX_N];
    copy_input(xr_rc, xi_rc, input, n);
    fft_1d_recursive(xr_rc, xi_rc, tr, ti, n);
    print_fft("Recursive FFT (all-ones input)", xr_rc, xi_rc, n);

    printf("  Expected: X[0]=%.1f, all others ~0\n", (float)n);

    /* Impulse signal: x[0]=1, rest=0 → all X[k]=1 */
    for (int i = 0; i < n; i++) input[i] = 0.0f;
    input[0] = 1.0f;

    generate_twiddle_factors(tr, ti, n);
    copy_input(xr_it, xi_it, input, n);
    fft_1d_iterative(xr_it, xi_it, tr, ti, n);
    print_fft("Iterative FFT (impulse input)", xr_it, xi_it, n);
    printf("  Expected: all X[k] = 1.0000 + 0.0000i\n");
}

/* ---------- Performance benchmarks ---------- */
typedef struct {
    int    n;
    unsigned long cycles_iter;
    unsigned long cycles_rec;
} PerfResult;

static PerfResult bench(int n) {
    PerfResult r;
    r.n = n;

    float input[MAX_N];
    for (int i = 0; i < n; i++) input[i] = (float)i;

    float tr[MAX_N/2], ti[MAX_N/2];
    float xr[MAX_N], xi[MAX_N];

    /* --- Iterative --- */
    generate_twiddle_factors(tr, ti, n);
    copy_input(xr, xi, input, n);
    unsigned long c0 = read_cycles();
    fft_1d_iterative(xr, xi, tr, ti, n);
    unsigned long c1 = read_cycles();
    r.cycles_iter = c1 - c0;

    /* --- Recursive --- */
    generate_twiddle_factors(tr, ti, n);
    copy_input(xr, xi, input, n);
    unsigned long c2 = read_cycles();
    fft_1d_recursive(xr, xi, tr, ti, n);
    unsigned long c3 = read_cycles();
    r.cycles_rec = c3 - c2;

    return r;
}

/* Stack memory estimates (bytes):
 *   butterfly_iterative frame: 80 bytes, called once
 *   fft_1d_iterative frame:    48 bytes
 *   bit_reverse_array frame:   48 bytes
 *   reverse_bits: no frame (leaf, uses t-regs only)
 *   butterfly_recursive: 64 bytes * depth = 64 * log2(n)
 *   fft_1d_recursive frame:    48 bytes
 *
 *   Iterative total stack = 80 + 48 + 48 = 176 bytes
 *   Recursive total stack = 48 + 48 + 64*log2(n) bytes
 *   Plus twiddle array on heap: n/2 * 2 * 4 bytes
 */
static int stack_iterative(int n) { return 176; }
static int stack_recursive(int n) {
    int log_n = 0; int tmp = n;
    while (tmp > 1) { tmp >>= 1; log_n++; }
    return 48 + 48 + 64 * log_n;
}

int main(void) {
    printf("============================================================\n");
    printf("  Milestone 2: 1D FFT in RISC-V Assembly\n");
    printf("============================================================\n");

    test_math();
    test_bit_reversal();
    test_array_reorder();

    printf("\n=== FFT Correctness (n=8) ===");
    test_fft(8);

    /* Performance table */
    int sizes[] = {8, 16, 32};
    PerfResult results[3];
    printf("\n\n=== Performance Benchmarks ===\n");
    for (int i = 0; i < 3; i++) {
        results[i] = bench(sizes[i]);
    }

    printf("\n+------+----------------+------------------+------------------+------------------+\n");
    printf("| N    | Iter Cycles    | Rec  Cycles      | Iter Stack(B)    | Rec  Stack(B)    |\n");
    printf("+------+----------------+------------------+------------------+------------------+\n");
    for (int i = 0; i < 3; i++) {
        int n = results[i].n;
        printf("| %-4d | %-14lu | %-16lu | %-16d | %-16d |\n",
               n,
               results[i].cycles_iter,
               results[i].cycles_rec,
               stack_iterative(n),
               stack_recursive(n));
    }
    printf("+------+----------------+------------------+------------------+------------------+\n");

    printf("\n=== Energy Analysis ===\n");
    printf("Fewer cycles = less CPU active time = less energy consumed.\n");
    printf("Less stack = lower memory power consumption.\n\n");
    for (int i = 0; i < 3; i++) {
        int n = results[i].n;
        const char *faster = (results[i].cycles_iter <= results[i].cycles_rec)
                             ? "Iterative" : "Recursive";
        const char *smaller = (stack_iterative(n) <= stack_recursive(n))
                             ? "Iterative" : "Recursive";
        printf("  n=%-2d: Fewer cycles -> %s | Less memory -> %s\n",
               n, faster, smaller);
    }
    printf("\nConclusion: Iterative FFT is preferred for battery-powered devices.\n");
    printf("  - Fixed stack depth avoids recursive call overhead.\n");
    printf("  - No repeated function-call prologue/epilogue overhead.\n");
    printf("  - Predictable memory access pattern → better cache utilization.\n");

    return 0;
}
