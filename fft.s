# =============================================================================
# fft.s  –  1D FFT in RISC-V (RV64GC) assembly
# Functions exported (callable from C):
#   void      my_sin(float x)              -> fa0
#   void      my_cos(float x)              -> fa0
#   void      generate_twiddle_factors(float *tr, float *ti, int n)
#   unsigned  reverse_bits(unsigned x, int log_n)  -> a0
#   void      bit_reverse_array(float *xr, float *xi, int n)
#   void      butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n)
#   void      fft_1d_iterative(float *xr, float *xi, float *tr, float *ti, int n)
#   void      butterfly_recursive(float *xr, float *xi, float *tr, float *ti, int n, int total_n)
#   void      fft_1d_recursive(float *xr, float *xi, float *tr, float *ti, int n)
# =============================================================================

    .section .rodata
PI:         .float  3.14159265358979323846
TWO_PI:     .float  6.28318530717958647692
NEG_TWO_PI: .float -6.28318530717958647692
ONE:        .float  1.0
NEG_ONE:    .float -1.0
HALF_PI:    .float  1.5707963267948966

    .section .text

# =============================================================================
# float my_sin(float x)
#   Taylor series: sin(x) = x - x^3/3! + x^5/5! - x^7/7! + ...  (10 terms)
#   Argument reduction to [-PI, PI] first.
# Calling convention: argument in fa0, return in fa0
# =============================================================================
    .globl my_sin
    .type  my_sin, @function
my_sin:
    addi    sp, sp, -32
    sd      ra, 24(sp)
    fsd     fs0, 16(sp)
    fsd     fs1,  8(sp)
    fsd     fs2,  0(sp)

    fmv.s   fs0, fa0                # fs0 = x (working copy)

    # --- Argument reduction: wrap to [-PI, PI] ---
    la      t0, TWO_PI
    flw     ft1, 0(t0)
    la      t0, PI
    flw     ft2, 0(t0)              # ft2 = PI

sin_reduce_pos:
    flt.s   t0, ft2, fs0            # if PI < x
    beqz    t0, sin_reduce_neg
    fsub.s  fs0, fs0, ft1           # x -= 2*PI
    j       sin_reduce_pos

sin_reduce_neg:
    fneg.s  ft3, ft2                # ft3 = -PI
    flt.s   t0, fs0, ft3            # if x < -PI
    beqz    t0, sin_taylor
    fadd.s  fs0, fs0, ft1           # x += 2*PI
    j       sin_reduce_neg

sin_taylor:
    # res = 0,  term = x,  x_sq = x*x
    fmv.w.x fs1, zero               # fs1 = res = 0.0
    fmv.s   fs2, fs0                # fs2 = term = x
    fmul.s  ft4, fs0, fs0           # ft4 = x_sq

    # i = 1 .. 10
    li      t1, 1                   # t1 = i
sin_loop:
    li      t2, 11
    bge     t1, t2, sin_done

    fadd.s  fs1, fs1, fs2           # res += term

    # next term: term *= -x_sq / ((2i)*(2i+1))
    # compute denominator as float
    slli    t3, t1, 1               # t3 = 2i
    addi    t4, t3, 1               # t4 = 2i+1
    mul     t5, t3, t4              # t5 = (2i)*(2i+1)
    fcvt.s.w ft5, t5                # ft5 = float(denom)
    fdiv.s  ft6, ft4, ft5           # ft6 = x_sq / denom
    fneg.s  ft6, ft6                # ft6 = -x_sq / denom
    fmul.s  fs2, fs2, ft6           # term *= (-x_sq/denom)

    addi    t1, t1, 1
    j       sin_loop

sin_done:
    fmv.s   fa0, fs1                # return res

    fld     fs2,  0(sp)
    fld     fs1,  8(sp)
    fld     fs0, 16(sp)
    ld      ra,  24(sp)
    addi    sp, sp, 32
    ret

# =============================================================================
# float my_cos(float x)
#   cos(x) = sin(x + PI/2)
# =============================================================================
    .globl my_cos
    .type  my_cos, @function
my_cos:
    la      t0, HALF_PI
    flw     ft0, 0(t0)
    fadd.s  fa0, fa0, ft0           # x += PI/2
    j       my_sin                  # tail-call my_sin

# =============================================================================
# void generate_twiddle_factors(float *twiddle_real, float *twiddle_imag, int n)
#   W_k = exp(-j * 2*PI*k / n)  for k = 0 .. n/2-1
#   real part = cos(-2*PI*k/n),  imag part = sin(-2*PI*k/n)
# a0 = twiddle_real,  a1 = twiddle_imag,  a2 = n
# =============================================================================
    .globl generate_twiddle_factors
    .type  generate_twiddle_factors, @function
generate_twiddle_factors:
    addi    sp, sp, -64
    sd      ra,  56(sp)
    sd      s0,  48(sp)
    sd      s1,  40(sp)
    sd      s2,  32(sp)
    sd      s3,  24(sp)
    fsd     fs0, 16(sp)
    fsd     fs1,  8(sp)

    mv      s0, a0                  # s0 = twiddle_real ptr
    mv      s1, a1                  # s1 = twiddle_imag ptr
    mv      s2, a2                  # s2 = n
    li      s3, 0                   # s3 = k = 0

    srli    t0, s2, 1               # half = n/2
    fcvt.s.w ft2, s2                # ft2 = float(n)
    la      t1, NEG_TWO_PI
    flw     fs0, 0(t1)              # fs0 = -2*PI

twiddle_loop:
    bge     s3, t0, twiddle_done    # k >= n/2 → stop

    fcvt.s.w fa0, s3                # fa0 = float(k)
    fmul.s  fa0, fs0, fa0           # fa0 = -2*PI*k
    fdiv.s  fa0, fa0, ft2           # fa0 = -2*PI*k/n  (= angle)

    # cos(angle) → twiddle_real[k]
    fsd     fa0, 0(sp)              # spill angle
    call    my_cos
    flw     ft3, 0(sp)              # restore angle (stored as float)
    # (my_cos trashed fa0 with result)
    fsw     fa0, 0(s0)              # twiddle_real[k] = cos(angle)
    addi    s0, s0, 4

    # sin(angle) → twiddle_imag[k]
    flw     fa0, 0(sp)              # reload angle
    call    my_sin
    fsw     fa0, 0(s1)              # twiddle_imag[k] = sin(angle)
    addi    s1, s1, 4

    addi    s3, s3, 1               # k++
    srli    t0, s2, 1               # recompute n/2 (register reuse)
    fcvt.s.w ft2, s2
    j       twiddle_loop

twiddle_done:
    fld     fs1,  8(sp)
    fld     fs0, 16(sp)
    ld      s3,  24(sp)
    ld      s2,  32(sp)
    ld      s1,  40(sp)
    ld      s0,  48(sp)
    ld      ra,  56(sp)
    addi    sp, sp, 64
    ret

# =============================================================================
# unsigned int reverse_bits(unsigned int x, int log_n)
#   a0 = x,  a1 = log_n
#   returns reversed value in a0
# =============================================================================
    .globl reverse_bits
    .type  reverse_bits, @function
reverse_bits:
    li      t0, 0                   # reversed = 0
    mv      t1, a1                  # loop counter = log_n

rb_loop:
    beqz    t1, rb_done
    slli    t0, t0, 1               # reversed <<= 1
    andi    t2, a0, 1               # t2 = x & 1
    or      t0, t0, t2              # reversed |= LSB
    srli    a0, a0, 1               # x >>= 1
    addi    t1, t1, -1
    j       rb_loop

rb_done:
    mv      a0, t0
    ret

# =============================================================================
# void bit_reverse_array(float *x_real, float *x_imag, int n)
#   a0 = x_real,  a1 = x_imag,  a2 = n
# =============================================================================
    .globl bit_reverse_array
    .type  bit_reverse_array, @function
bit_reverse_array:
    addi    sp, sp, -48
    sd      ra,  40(sp)
    sd      s0,  32(sp)
    sd      s1,  24(sp)
    sd      s2,  16(sp)
    sd      s3,   8(sp)
    sd      s4,   0(sp)

    mv      s0, a0                  # s0 = x_real
    mv      s1, a1                  # s1 = x_imag
    mv      s2, a2                  # s2 = n
    li      s3, 0                   # s3 = i = 0

    # compute log_n
    mv      a0, s2
    call    log2_int                # returns in a0
    mv      s4, a0                  # s4 = log_n

bra_loop:
    bge     s3, s2, bra_done        # i >= n → stop

    mv      a0, s3
    mv      a1, s4
    call    reverse_bits            # j = reverse_bits(i, log_n)
    mv      t0, a0                  # t0 = j

    ble     t0, s3, bra_next        # only swap if i < j

    # swap x_real[i] and x_real[j]
    slli    t1, s3, 2
    add     t1, s0, t1              # &x_real[i]
    slli    t2, t0, 2
    add     t2, s0, t2              # &x_real[j]
    flw     ft0, 0(t1)
    flw     ft1, 0(t2)
    fsw     ft1, 0(t1)
    fsw     ft0, 0(t2)

    # swap x_imag[i] and x_imag[j]
    slli    t1, s3, 2
    add     t1, s1, t1
    slli    t2, t0, 2
    add     t2, s1, t2
    flw     ft0, 0(t1)
    flw     ft1, 0(t2)
    fsw     ft1, 0(t1)
    fsw     ft0, 0(t2)

bra_next:
    addi    s3, s3, 1
    j       bra_loop

bra_done:
    ld      s4,   0(sp)
    ld      s3,   8(sp)
    ld      s2,  16(sp)
    ld      s1,  24(sp)
    ld      s0,  32(sp)
    ld      ra,  40(sp)
    addi    sp, sp, 48
    ret

# =============================================================================
# int log2_int(int n)   -- helper used by bit_reverse_array
# =============================================================================
    .globl log2_int
    .type  log2_int, @function
log2_int:
    li      a1, 0
l2_loop:
    li      t0, 1
    ble     a0, t0, l2_done
    srli    a0, a0, 1
    addi    a1, a1, 1
    j       l2_loop
l2_done:
    mv      a0, a1
    ret

# =============================================================================
# void butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n)
#   a0=xr, a1=xi, a2=tr, a3=ti, a4=n
# =============================================================================
    .globl butterfly_iterative
    .type  butterfly_iterative, @function
butterfly_iterative:
    addi    sp, sp, -80
    sd      ra,  72(sp)
    sd      s0,  64(sp)
    sd      s1,  56(sp)
    sd      s2,  48(sp)
    sd      s3,  40(sp)
    sd      s4,  32(sp)
    sd      s5,  24(sp)
    sd      s6,  16(sp)
    sd      s7,   8(sp)
    sd      s8,   0(sp)

    mv      s0, a0                  # xr
    mv      s1, a1                  # xi
    mv      s2, a2                  # tr
    mv      s3, a3                  # ti
    mv      s4, a4                  # n

    # stages = log2(n)
    mv      a0, s4
    call    log2_int
    mv      s5, a0                  # s5 = stages

    li      s6, 1                   # s6 = stage = 1

bi_stage_loop:
    bgt     s6, s5, bi_done         # stage > stages → done

    li      t0, 1
    sll     t0, t0, s6              # t0 = m = 1 << stage
    srli    t1, t0, 1               # t1 = half_m = m/2
    div     t2, s4, t0              # t2 = twiddle_step = n/m

    li      s7, 0                   # s7 = k = 0

bi_k_loop:
    bge     s7, s4, bi_stage_next   # k >= n → next stage

    li      s8, 0                   # s8 = j = 0

bi_j_loop:
    bge     s8, t1, bi_k_next       # j >= half_m → next k

    # top = k+j,  bottom = k+j+half_m
    add     t3, s7, s8              # t3 = top index
    add     t4, t3, t1              # t4 = bottom index

    # twiddle index = j * twiddle_step
    mul     t5, s8, t2              # t5 = twiddle index

    # load twiddle: wr = tr[t5], wi = ti[t5]
    slli    t6, t5, 2
    add     t6, s2, t6
    flw     ft0, 0(t6)              # ft0 = wr

    slli    t6, t5, 2
    add     t6, s3, t6
    flw     ft1, 0(t6)              # ft1 = wi

    # load bottom element
    slli    t6, t4, 2
    add     t6, s0, t6
    flw     ft2, 0(t6)              # ft2 = xr[bottom]

    slli    t6, t4, 2
    add     t6, s1, t6
    flw     ft3, 0(t6)              # ft3 = xi[bottom]

    # tr_val = wr*xr[bot] - wi*xi[bot]
    fmul.s  ft4, ft0, ft2
    fmul.s  ft5, ft1, ft3
    fsub.s  ft4, ft4, ft5           # ft4 = tr_val

    # ti_val = wr*xi[bot] + wi*xr[bot]
    fmul.s  ft5, ft0, ft3
    fmul.s  ft6, ft1, ft2
    fadd.s  ft5, ft5, ft6           # ft5 = ti_val

    # load top element
    slli    t6, t3, 2
    add     t6, s0, t6
    flw     ft6, 0(t6)              # ft6 = xr[top]
    slli    t6, t3, 2
    add     t6, s1, t6
    flw     ft7, 0(t6)              # ft7 = xi[top]

    # xr[bot] = xr[top] - tr_val
    fsub.s  fa0, ft6, ft4
    slli    t6, t4, 2
    add     t6, s0, t6
    fsw     fa0, 0(t6)

    # xi[bot] = xi[top] - ti_val
    fsub.s  fa0, ft7, ft5
    slli    t6, t4, 2
    add     t6, s1, t6
    fsw     fa0, 0(t6)

    # xr[top] = xr[top] + tr_val
    fadd.s  fa0, ft6, ft4
    slli    t6, t3, 2
    add     t6, s0, t6
    fsw     fa0, 0(t6)

    # xi[top] = xi[top] + ti_val
    fadd.s  fa0, ft7, ft5
    slli    t6, t3, 2
    add     t6, s1, t6
    fsw     fa0, 0(t6)

    addi    s8, s8, 1               # j++
    j       bi_j_loop

bi_k_next:
    add     s7, s7, t0              # k += m
    j       bi_k_loop

bi_stage_next:
    addi    s6, s6, 1               # stage++
    # recompute m, half_m, twiddle_step for next stage
    li      t0, 1
    sll     t0, t0, s6
    srli    t1, t0, 1
    div     t2, s4, t0
    j       bi_stage_loop

bi_done:
    ld      s8,   0(sp)
    ld      s7,   8(sp)
    ld      s6,  16(sp)
    ld      s5,  24(sp)
    ld      s4,  32(sp)
    ld      s3,  40(sp)
    ld      s2,  48(sp)
    ld      s1,  56(sp)
    ld      s0,  64(sp)
    ld      ra,  72(sp)
    addi    sp, sp, 80
    ret

# =============================================================================
# void fft_1d_iterative(float *xr, float *xi, float *tr, float *ti, int n)
# =============================================================================
    .globl fft_1d_iterative
    .type  fft_1d_iterative, @function
fft_1d_iterative:
    addi    sp, sp, -48
    sd      ra,  40(sp)
    sd      s0,  32(sp)
    sd      s1,  24(sp)
    sd      s2,  16(sp)
    sd      s3,   8(sp)
    sd      s4,   0(sp)

    mv      s0, a0
    mv      s1, a1
    mv      s2, a2
    mv      s3, a3
    mv      s4, a4

    # bit_reverse_array(xr, xi, n)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s4
    call    bit_reverse_array

    # butterfly_iterative(xr, xi, tr, ti, n)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s2
    mv      a3, s3
    mv      a4, s4
    call    butterfly_iterative

    ld      s4,   0(sp)
    ld      s3,   8(sp)
    ld      s2,  16(sp)
    ld      s1,  24(sp)
    ld      s0,  32(sp)
    ld      ra,  40(sp)
    addi    sp, sp, 48
    ret

# =============================================================================
# void butterfly_recursive(float *xr, float *xi, float *tr, float *ti, int n, int total_n)
#   a0=xr, a1=xi, a2=tr, a3=ti, a4=n, a5=total_n
# =============================================================================
    .globl butterfly_recursive
    .type  butterfly_recursive, @function
butterfly_recursive:
    # base case: n <= 1, return
    li      t0, 1
    ble     a4, t0, br_base

    addi    sp, sp, -64
    sd      ra,  56(sp)
    sd      s0,  48(sp)
    sd      s1,  40(sp)
    sd      s2,  32(sp)
    sd      s3,  24(sp)
    sd      s4,  16(sp)
    sd      s5,   8(sp)
    sd      s6,   0(sp)

    mv      s0, a0                  # xr
    mv      s1, a1                  # xi
    mv      s2, a2                  # tr
    mv      s3, a3                  # ti
    mv      s4, a4                  # n
    mv      s5, a5                  # total_n
    srli    s6, s4, 1               # s6 = half = n/2

    # recursive call on even (first half): butterfly_recursive(xr, xi, tr, ti, half, total_n)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s2
    mv      a3, s3
    mv      a4, s6
    mv      a5, s5
    call    butterfly_recursive

    # recursive call on odd (second half): butterfly_recursive(xr+half, xi+half, tr, ti, half, total_n)
    slli    t0, s6, 2
    add     a0, s0, t0              # xr + half*4
    add     a1, s1, t0              # xi + half*4
    mv      a2, s2
    mv      a3, s3
    mv      a4, s6
    mv      a5, s5
    call    butterfly_recursive

    # combine: for j = 0 .. half-1
    #   t_idx = j * (total_n / n)
    #   wr = tr[t_idx], wi = ti[t_idx]
    #   tr_val = xr[j+half]*wr - xi[j+half]*wi
    #   ti_val = xr[j+half]*wi + xi[j+half]*wr
    #   xr[j+half] = xr[j] - tr_val
    #   xi[j+half] = xi[j] - ti_val
    #   xr[j]      = xr[j] + tr_val
    #   xi[j]      = xi[j] + ti_val

    div     t2, s5, s4              # t2 = total_n / n
    li      t3, 0                   # t3 = j = 0

br_combine:
    bge     t3, s6, br_combine_done

    mul     t4, t3, t2              # t4 = j * (total_n/n) = twiddle index

    slli    t5, t4, 2
    add     t5, s2, t5
    flw     ft0, 0(t5)              # ft0 = wr

    slli    t5, t4, 2
    add     t5, s3, t5
    flw     ft1, 0(t5)              # ft1 = wi

    # bottom = j + half
    add     t5, t3, s6

    slli    t6, t5, 2
    add     t6, s0, t6
    flw     ft2, 0(t6)              # ft2 = xr[j+half]

    slli    t6, t5, 2
    add     t6, s1, t6
    flw     ft3, 0(t6)              # ft3 = xi[j+half]

    # tr_val = xr[j+half]*wr - xi[j+half]*wi
    fmul.s  ft4, ft2, ft0
    fmul.s  ft5, ft3, ft1
    fsub.s  ft4, ft4, ft5

    # ti_val = xr[j+half]*wi + xi[j+half]*wr
    fmul.s  ft5, ft2, ft1
    fmul.s  ft6, ft3, ft0
    fadd.s  ft5, ft5, ft6

    # load xr[j], xi[j]
    slli    t6, t3, 2
    add     t6, s0, t6
    flw     ft6, 0(t6)              # ft6 = xr[j]
    slli    t6, t3, 2
    add     t6, s1, t6
    flw     ft7, 0(t6)              # ft7 = xi[j]

    # xr[j+half] = xr[j] - tr_val
    fsub.s  fa5, ft6, ft4
    add     t6, t3, s6
    slli    t6, t6, 2
    add     t6, s0, t6
    fsw     fa5, 0(t6)

    # xi[j+half] = xi[j] - ti_val
    fsub.s  fa5, ft7, ft5
    add     t6, t3, s6
    slli    t6, t6, 2
    add     t6, s1, t6
    fsw     fa5, 0(t6)

    # xr[j] = xr[j] + tr_val
    fadd.s  fa5, ft6, ft4
    slli    t6, t3, 2
    add     t6, s0, t6
    fsw     fa5, 0(t6)

    # xi[j] = xi[j] + ti_val
    fadd.s  fa5, ft7, ft5
    slli    t6, t3, 2
    add     t6, s1, t6
    fsw     fa5, 0(t6)

    addi    t3, t3, 1               # j++
    j       br_combine

br_combine_done:
    ld      s6,   0(sp)
    ld      s5,   8(sp)
    ld      s4,  16(sp)
    ld      s3,  24(sp)
    ld      s2,  32(sp)
    ld      s1,  40(sp)
    ld      s0,  48(sp)
    ld      ra,  56(sp)
    addi    sp, sp, 64
    ret

br_base:
    ret

# =============================================================================
# void fft_1d_recursive(float *xr, float *xi, float *tr, float *ti, int n)
# =============================================================================
    .globl fft_1d_recursive
    .type  fft_1d_recursive, @function
fft_1d_recursive:
    addi    sp, sp, -48
    sd      ra,  40(sp)
    sd      s0,  32(sp)
    sd      s1,  24(sp)
    sd      s2,  16(sp)
    sd      s3,   8(sp)
    sd      s4,   0(sp)

    mv      s0, a0
    mv      s1, a1
    mv      s2, a2
    mv      s3, a3
    mv      s4, a4

    # bit_reverse_array(xr, xi, n)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s4
    call    bit_reverse_array

    # butterfly_recursive(xr, xi, tr, ti, n, n)
    mv      a0, s0
    mv      a1, s1
    mv      a2, s2
    mv      a3, s3
    mv      a4, s4
    mv      a5, s4                  # total_n = n
    call    butterfly_recursive

    ld      s4,   0(sp)
    ld      s3,   8(sp)
    ld      s2,  16(sp)
    ld      s1,  24(sp)
    ld      s0,  32(sp)
    ld      ra,  40(sp)
    addi    sp, sp, 48
    ret
