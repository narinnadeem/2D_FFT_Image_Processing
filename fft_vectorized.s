# =============================================================================
# fft_vectorized.s  -  1D FFT using RISC-V Vector (RVV) instructions
# Requires: -march=rv64gcv  (assembler)  and  -cpu rv64,v=true,vlen=128  (QEMU)
#
# Functions exported:
#   void vec_generate_twiddle(float *tr, float *ti, int n)
#   void vec_bit_reverse_array(float *xr, float *xi, int n)
#   void vec_butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n)
#   void fft_1d_vectorized(float *xr, float *xi, float *tr, float *ti, int n)
#
# Reused from fft.s (scalar helpers):
#   my_sin, my_cos, log2_int, reverse_bits
# =============================================================================

    .section .rodata
VEC_PI:         .float  3.14159265358979323846
VEC_TWO_PI:     .float  6.28318530717958647692
VEC_NEG_TWO_PI: .float -6.28318530717958647692
VEC_HALF_PI:    .float  1.5707963267948966

# Index array for twiddle factor generation: 0.0, 1.0, 2.0, ... 31.0
# The vector code loads chunks of this to compute angles in parallel.
    .section .data
    .align 4
vec_index_array:
    .float 0.0,  1.0,  2.0,  3.0,  4.0,  5.0,  6.0,  7.0
    .float 8.0,  9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0
    .float 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0
    .float 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0
    .float 32.0, 33.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0
    .float 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0
    .float 48.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 55.0
    .float 56.0, 57.0, 58.0, 59.0, 60.0, 61.0, 62.0, 63.0
    .float 64.0, 65.0, 66.0, 67.0, 68.0, 69.0, 70.0, 71.0
    .float 72.0, 73.0, 74.0, 75.0, 76.0, 77.0, 78.0, 79.0
    .float 80.0, 81.0, 82.0, 83.0, 84.0, 85.0, 86.0, 87.0
    .float 88.0, 89.0, 90.0, 91.0, 92.0, 93.0, 94.0, 95.0
    .float 96.0, 97.0, 98.0, 99.0,100.0,101.0,102.0,103.0
    .float 104.0,105.0,106.0,107.0,108.0,109.0,110.0,111.0
    .float 112.0,113.0,114.0,115.0,116.0,117.0,118.0,119.0
    .float 120.0,121.0,122.0,123.0,124.0,125.0,126.0,127.0
    .float 128.0,129.0,130.0,131.0,132.0,133.0,134.0,135.0
    .float 136.0,137.0,138.0,139.0,140.0,141.0,142.0,143.0
    .float 144.0,145.0,146.0,147.0,148.0,149.0,150.0,151.0
    .float 152.0,153.0,154.0,155.0,156.0,157.0,158.0,159.0
    .float 160.0,161.0,162.0,163.0,164.0,165.0,166.0,167.0
    .float 168.0,169.0,170.0,171.0,172.0,173.0,174.0,175.0
    .float 176.0,177.0,178.0,179.0,180.0,181.0,182.0,183.0
    .float 184.0,185.0,186.0,187.0,188.0,189.0,190.0,191.0
    .float 192.0,193.0,194.0,195.0,196.0,197.0,198.0,199.0
    .float 200.0,201.0,202.0,203.0,204.0,205.0,206.0,207.0
    .float 208.0,209.0,210.0,211.0,212.0,213.0,214.0,215.0
    .float 216.0,217.0,218.0,219.0,220.0,221.0,222.0,223.0
    .float 224.0,225.0,226.0,227.0,228.0,229.0,230.0,231.0
    .float 232.0,233.0,234.0,235.0,236.0,237.0,238.0,239.0
    .float 240.0,241.0,242.0,243.0,244.0,245.0,246.0,247.0
    .float 248.0,249.0,250.0,251.0,252.0,253.0,254.0,255.0
    .float 256.0,257.0,258.0,259.0,260.0,261.0,262.0,263.0
    .float 264.0,265.0,266.0,267.0,268.0,269.0,270.0,271.0
    .float 272.0,273.0,274.0,275.0,276.0,277.0,278.0,279.0
    .float 280.0,281.0,282.0,283.0,284.0,285.0,286.0,287.0
    .float 288.0,289.0,290.0,291.0,292.0,293.0,294.0,295.0
    .float 296.0,297.0,298.0,299.0,300.0,301.0,302.0,303.0
    .float 304.0,305.0,306.0,307.0,308.0,309.0,310.0,311.0
    .float 312.0,313.0,314.0,315.0,316.0,317.0,318.0,319.0
    .float 320.0,321.0,322.0,323.0,324.0,325.0,326.0,327.0
    .float 328.0,329.0,330.0,331.0,332.0,333.0,334.0,335.0
    .float 336.0,337.0,338.0,339.0,340.0,341.0,342.0,343.0
    .float 344.0,345.0,346.0,347.0,348.0,349.0,350.0,351.0
    .float 352.0,353.0,354.0,355.0,356.0,357.0,358.0,359.0
    .float 360.0,361.0,362.0,363.0,364.0,365.0,366.0,367.0
    .float 368.0,369.0,370.0,371.0,372.0,373.0,374.0,375.0
    .float 376.0,377.0,378.0,379.0,380.0,381.0,382.0,383.0
    .float 384.0,385.0,386.0,387.0,388.0,389.0,390.0,391.0
    .float 392.0,393.0,394.0,395.0,396.0,397.0,398.0,399.0
    .float 400.0,401.0,402.0,403.0,404.0,405.0,406.0,407.0
    .float 408.0,409.0,410.0,411.0,412.0,413.0,414.0,415.0
    .float 416.0,417.0,418.0,419.0,420.0,421.0,422.0,423.0
    .float 424.0,425.0,426.0,427.0,428.0,429.0,430.0,431.0
    .float 432.0,433.0,434.0,435.0,436.0,437.0,438.0,439.0
    .float 440.0,441.0,442.0,443.0,444.0,445.0,446.0,447.0
    .float 448.0,449.0,450.0,451.0,452.0,453.0,454.0,455.0
    .float 456.0,457.0,458.0,459.0,460.0,461.0,462.0,463.0
    .float 464.0,465.0,466.0,467.0,468.0,469.0,470.0,471.0
    .float 472.0,473.0,474.0,475.0,476.0,477.0,478.0,479.0
    .float 480.0,481.0,482.0,483.0,484.0,485.0,486.0,487.0
    .float 488.0,489.0,490.0,491.0,492.0,493.0,494.0,495.0
    .float 496.0,497.0,498.0,499.0,500.0,501.0,502.0,503.0
    .float 504.0,505.0,506.0,507.0,508.0,509.0,510.0,511.0

    .section .text

# =============================================================================
# void vec_generate_twiddle(float *tr, float *ti, int n)
#
# Generates twiddle factors W_k = e^(-j*2*PI*k/n) for k = 0 .. n/2-1
# in chunks using vector loads from the index array.
#
# Key vector instructions used:
#   vsetvli  - set vector length based on available registers and element type
#   vle32.v  - unit-stride vector load (32-bit floats)
#   vfmul.vf - vector-scalar float multiply
#   vfdiv.vf - vector-scalar float divide
#   vse32.v  - unit-stride vector store
#
# The sin/cos calls remain scalar because vectorizing transcendental functions
# requires lookup tables (beyond scope); we vectorize the index arithmetic
# and memory operations around them.
# =============================================================================
    .globl vec_generate_twiddle
    .type  vec_generate_twiddle, @function
vec_generate_twiddle:
    # Stack: 96 bytes
    # 88=ra,80=s0,72=s1,64=s2,56=s3,48=s4,40=s5,32=s6,24=s7,16=fs0,8=fs1,0=angle(float)
    addi    sp, sp, -96
    sd      ra,  88(sp)
    sd      s0,  80(sp)
    sd      s1,  72(sp)
    sd      s2,  64(sp)
    sd      s3,  56(sp)
    sd      s4,  48(sp)
    sd      s5,  40(sp)
    sd      s6,  32(sp)
    sd      s7,  24(sp)
    fsd     fs0, 16(sp)
    fsd     fs1,  8(sp)
    # 0(sp) = 4-byte angle spill slot

    mv      s0, a0                  # tr ptr
    mv      s1, a1                  # ti ptr
    mv      s2, a2                  # n
    srli    s3, s2, 1               # half = n/2
    li      s4, 0                   # k = 0

    la      t0, VEC_NEG_TWO_PI
    flw     fs0, 0(t0)              # fs0 = -2*PI (saved, survives calls)
    fcvt.s.w fs1, s2                # fs1 = float(n) (saved, survives calls)
    la      s7, vec_index_array     # s7 = index array base (saved)

vec_twiddle_chunk:
    bge     s4, s3, vec_twiddle_done

    sub     t0, s3, s4
    vsetvli s5, t0, e32, m1, ta, ma # s5 = vl (saved reg - survives calls!)

    slli    t0, s4, 2
    add     t0, s7, t0
    vle32.v v0, (t0)                # v0 = float indices [k..k+vl-1]

    vfmul.vf v1, v0, fs0            # v1 = -2*PI * index
    vfdiv.vf v1, v1, fs1            # v1 = angles

    li      s6, 0                   # s6 = elem = 0 (saved reg - survives calls!)

vec_twiddle_elem:
    bge     s6, s5, vec_twiddle_chunk_done

    vslidedown.vx v2, v1, s6        # v2[0] = angle[s6]
    vfmv.f.s fa0, v2                # fa0 = angle
    fsw     fa0, 0(sp)              # spill angle (4 bytes at 0(sp))

    call    my_cos                  # fa0 = cos(angle); s regs all preserved by callee

    add     t0, s4, s6              # index = k + elem
    slli    t0, t0, 2
    add     t0, s0, t0
    fsw     fa0, 0(t0)              # tr[k+elem] = cos

    flw     fa0, 0(sp)              # reload angle
    call    my_sin                  # fa0 = sin(angle)

    add     t0, s4, s6
    slli    t0, t0, 2
    add     t0, s1, t0
    fsw     fa0, 0(t0)              # ti[k+elem] = sin

    addi    s6, s6, 1
    j       vec_twiddle_elem

vec_twiddle_chunk_done:
    add     s4, s4, s5              # k += vl
    j       vec_twiddle_chunk

vec_twiddle_done:
    fld     fs1,  8(sp)
    fld     fs0, 16(sp)
    ld      s7,  24(sp)
    ld      s6,  32(sp)
    ld      s5,  40(sp)
    ld      s4,  48(sp)
    ld      s3,  56(sp)
    ld      s2,  64(sp)
    ld      s1,  72(sp)
    ld      s0,  80(sp)
    ld      ra,  88(sp)
    addi    sp, sp, 96
    ret

# =============================================================================
# void vec_bit_reverse_array(float *xr, float *xi, int n)
#
# Computes bit-reversed indices for the entire array using vector instructions,
# then reorders xr and xi using gather loads into a temporary buffer.
#
# Key vector instructions:
#   vid.v    - generate vector index sequence [0, 1, 2, ..., vl-1]
#   vsrl.vx  - vector shift right logical by scalar
#   vsll.vx  - vector shift left logical by scalar
#   vand.vx  - vector AND with scalar
#   vor.vv   - vector OR
#   vluxei32.v - indexed (gather) load: load xr[index_vec[i]] for each i
#   vse32.v  - store result back
# =============================================================================
    .globl vec_bit_reverse_array
    .type  vec_bit_reverse_array, @function
vec_bit_reverse_array:
    addi    sp, sp, -80
    sd      ra,  72(sp)
    sd      s0,  64(sp)
    sd      s1,  56(sp)
    sd      s2,  48(sp)
    sd      s3,  40(sp)
    sd      s4,  32(sp)
    sd      s5,  24(sp)

    mv      s0, a0                  # s0 = xr
    mv      s1, a1                  # s1 = xi
    mv      s2, a2                  # s2 = n

    # Compute log_n
    mv      a0, s2
    call    log2_int
    mv      s3, a0                  # s3 = log_n

    # Allocate temp buffers on heap via stack pointer trick
    # We need n*4 bytes for each of: rev_idx (int), tmp_xr, tmp_xi
    # Use alloca-style: adjust sp
    slli    t0, s2, 2               # t0 = n * 4 bytes
    sub     sp, sp, t0              # tmp_xi buffer
    mv      s5, sp
    sub     sp, sp, t0              # tmp_xr buffer
    mv      s4, sp
    sub     sp, sp, t0              # rev_idx buffer (int32)
    mv      t2, sp                  # t2 = rev_idx base

    # --- Step 1: Compute bit-reversed indices using vector instructions ---
    # Process in chunks of vl
    li      t3, 0                   # t3 = i = 0

vec_bitrev_idx_loop:
    bge     t3, s2, vec_bitrev_idx_done

    sub     t4, s2, t3              # remaining
    vsetvli t0, t4, e32, m1, ta, ma # set vl

    # Generate index sequence: v0 = [i, i+1, ..., i+vl-1]
    vid.v   v0                      # v0 = [0, 1, ..., vl-1]
    vmv.v.x v1, t3                  # v1 = [i, i, ..., i]
    vadd.vv v0, v0, v1              # v0 = [i, i+1, ..., i+vl-1]

    # Bit-reverse each index using scalar loop log_n times
    # reversed = 0; for b in range(log_n): reversed = (reversed<<1)|(x&1); x>>=1
    vmv.v.x v2, zero                # v2 = reversed = [0, 0, ...]
    vmv.v.v v3, v0                  # v3 = x = copy of indices
    mv      t5, s3                  # t5 = bit counter = log_n

vec_bitrev_bits:
    beqz    t5, vec_bitrev_bits_done
    li      t4, 1
    vsll.vx v2, v2, t4              # reversed <<= 1
    vand.vx v4, v3, t4              # v4 = x & 1
    vor.vv  v2, v2, v4              # reversed |= LSB
    vsrl.vx v3, v3, t4              # x >>= 1
    addi    t5, t5, -1
    j       vec_bitrev_bits

vec_bitrev_bits_done:
    # Store reversed indices to rev_idx buffer
    slli    t6, t3, 2
    add     t6, t2, t6              # &rev_idx[i]
    vse32.v v2, (t6)

    add     t3, t3, t0              # i += vl
    j       vec_bitrev_idx_loop

vec_bitrev_idx_done:

    # --- Step 2: Gather load xr and xi using reversed indices ---
    # tmp_xr[i] = xr[rev_idx[i]],  tmp_xi[i] = xi[rev_idx[i]]
    li      t3, 0                   # i = 0

vec_gather_loop:
    bge     t3, s2, vec_gather_done

    sub     t4, s2, t3
    vsetvli t0, t4, e32, m1, ta, ma

    # Load rev_idx[i..i+vl-1] into v2 (byte offsets = rev_idx * 4)
    slli    t6, t3, 2
    add     t6, t2, t6
    vle32.v v2, (t6)                # v2 = rev_idx[i..i+vl-1]

    # Convert indices to byte offsets (* 4)
    li      t4, 2
    vsll.vx v3, v2, t4              # v3 = rev_idx * 4 (byte offsets)

    # Gather load: tmp_xr[i] = xr[rev_idx[i]]
    vluxei32.v v4, (s0), v3         # v4 = xr[rev_idx[i..i+vl-1]]
    slli    t6, t3, 2
    add     t6, s4, t6
    vse32.v v4, (t6)                # tmp_xr[i..i+vl-1] = v4

    # Gather load: tmp_xi[i] = xi[rev_idx[i]]
    vluxei32.v v5, (s1), v3         # v5 = xi[rev_idx[i..i+vl-1]]
    slli    t6, t3, 2
    add     t6, s5, t6
    vse32.v v5, (t6)                # tmp_xi[i..i+vl-1] = v5

    add     t3, t3, t0
    j       vec_gather_loop

vec_gather_done:
    # --- Step 3: Copy tmp buffers back to xr, xi ---
    li      t3, 0

vec_copy_back:
    bge     t3, s2, vec_copy_done

    sub     t4, s2, t3
    vsetvli t0, t4, e32, m1, ta, ma

    slli    t6, t3, 2
    add     t6, s4, t6
    vle32.v v0, (t6)                # load tmp_xr[i..i+vl-1]
    slli    t6, t3, 2
    add     t6, s0, t6
    vse32.v v0, (t6)                # store to xr[i..i+vl-1]

    slli    t6, t3, 2
    add     t6, s5, t6
    vle32.v v1, (t6)                # load tmp_xi[i..i+vl-1]
    slli    t6, t3, 2
    add     t6, s1, t6
    vse32.v v1, (t6)                # store to xi[i..i+vl-1]

    add     t3, t3, t0
    j       vec_copy_back

vec_copy_done:
    # Restore dynamic stack allocation
    slli    t0, s2, 2
    add     sp, sp, t0              # free rev_idx
    add     sp, sp, t0              # free tmp_xr
    add     sp, sp, t0              # free tmp_xi

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
# void vec_butterfly_iterative(float *xr, float *xi, float *tr, float *ti, int n)
#
# Vectorized iterative Cooley-Tukey butterfly.
# In each stage, the inner j-loop (which iterates over pairs within a group)
# is vectorized: we load half_m twiddle factors and half_m bottom elements
# at once, perform the butterfly computation on entire vectors, and store back.
#
# For large N where half_m > vl, we loop through the butterfly pairs in chunks.
#
# Key vector instructions:
#   vle32.v    - load twiddle factors (contiguous)
#   vlse32.v   - strided load for twiddle factors (stride = twiddle_step * 4)
#   vfmul.vv   - vector-vector float multiply
#   vfadd.vv   - vector-vector float add
#   vfsub.vv   - vector-vector float subtract
#   vse32.v    - store results
# =============================================================================
    .globl vec_butterfly_iterative
    .type  vec_butterfly_iterative, @function
vec_butterfly_iterative:
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

    mv      a0, s4
    call    log2_int
    mv      s5, a0                  # s5 = stages = log2(n)

    li      s6, 1                   # s6 = stage = 1

vbi_stage_loop:
    bgt     s6, s5, vbi_done

    li      t0, 1
    sll     t0, t0, s6              # t0 = m = 1 << stage
    srli    t1, t0, 1               # t1 = half_m = m / 2
    div     t2, s4, t0              # t2 = twiddle_step = n / m

    li      s7, 0                   # s7 = k = 0  (group start)

vbi_k_loop:
    bge     s7, s4, vbi_stage_next

    # Process j = 0 .. half_m-1 in vector chunks
    li      s8, 0                   # s8 = j_chunk_start = 0

vbi_j_chunk:
    bge     s8, t1, vbi_k_next      # j_chunk_start >= half_m -> next k

    # Remaining pairs in this group
    sub     t3, t1, s8              # remaining = half_m - j_chunk_start
    vsetvli t6, t3, e32, m1, ta, ma # t6 = vl

    # --- Load twiddle factors using strided load ---
    # twiddle index for j = j_chunk_start + [0..vl-1]  is
    #   (j_chunk_start + i) * twiddle_step
    # stride in bytes = twiddle_step * 4
    mul     t3, s8, t2              # base twiddle idx = j_chunk_start * twiddle_step
    slli    t3, t3, 2               # byte offset
    add     t3, s2, t3              # &tr[base_twiddle_idx]
    slli    t4, t2, 2               # stride = twiddle_step * 4 bytes
    vlse32.v v0, (t3), t4           # v0 = wr[j_chunk..j_chunk+vl-1] (strided)

    add     t3, s3, t3
    sub     t3, t3, s2              # same offset but in ti array
    mul     t3, s8, t2
    slli    t3, t3, 2
    add     t3, s3, t3
    vlse32.v v1, (t3), t4           # v1 = wi[j_chunk..j_chunk+vl-1] (strided)

    # --- Load bottom elements xr[k + j_chunk .. k + j_chunk + vl - 1 + half_m] ---
    add     t3, s7, s8              # base bottom = k + j_chunk_start
    add     t3, t3, t1              # + half_m
    slli    t3, t3, 2
    add     t3, s0, t3              # &xr[bottom_base]
    vle32.v v2, (t3)                # v2 = xr_bot[0..vl-1]

    add     t3, s7, s8
    add     t3, t3, t1
    slli    t3, t3, 2
    add     t3, s1, t3              # &xi[bottom_base]
    vle32.v v3, (t3)                # v3 = xi_bot[0..vl-1]

    # --- Load top elements xr[k + j_chunk .. k + j_chunk + vl - 1] ---
    add     t3, s7, s8              # base top = k + j_chunk_start
    slli    t3, t3, 2
    add     t3, s0, t3              # &xr[top_base]
    vle32.v v6, (t3)                # v6 = xr_top[0..vl-1]

    add     t3, s7, s8
    slli    t3, t3, 2
    add     t3, s1, t3              # &xi[top_base]
    vle32.v v7, (t3)                # v7 = xi_top[0..vl-1]

    # --- Butterfly computation (vectorized) ---
    # tr_val = wr * xr_bot - wi * xi_bot
    vfmul.vv v4, v0, v2             # v4 = wr * xr_bot
    vfmul.vv v5, v1, v3             # v5 = wi * xi_bot
    vfsub.vv v4, v4, v5             # v4 = tr_val

    # ti_val = wr * xi_bot + wi * xr_bot
    vfmul.vv v5, v0, v3             # v5 = wr * xi_bot
    vfmul.vv v8, v1, v2             # v8 = wi * xr_bot
    vfadd.vv v5, v5, v8             # v5 = ti_val

    # xr_top_new = xr_top + tr_val
    vfadd.vv v8, v6, v4
    # xr_bot_new = xr_top - tr_val
    vfsub.vv v9, v6, v4

    # xi_top_new = xi_top + ti_val
    vfadd.vv v10, v7, v5
    # xi_bot_new = xi_top - ti_val
    vfsub.vv v11, v7, v5

    # --- Store results back ---
    # top elements
    add     t3, s7, s8
    slli    t3, t3, 2
    add     t3, s0, t3
    vse32.v v8, (t3)                # xr[top] = xr_top_new

    add     t3, s7, s8
    slli    t3, t3, 2
    add     t3, s1, t3
    vse32.v v10, (t3)               # xi[top] = xi_top_new

    # bottom elements
    add     t3, s7, s8
    add     t3, t3, t1
    slli    t3, t3, 2
    add     t3, s0, t3
    vse32.v v9, (t3)                # xr[bot] = xr_bot_new

    add     t3, s7, s8
    add     t3, t3, t1
    slli    t3, t3, 2
    add     t3, s1, t3
    vse32.v v11, (t3)               # xi[bot] = xi_bot_new

    add     s8, s8, t6              # j_chunk_start += vl
    j       vbi_j_chunk

vbi_k_next:
    add     s7, s7, t0              # k += m
    j       vbi_k_loop

vbi_stage_next:
    addi    s6, s6, 1               # stage++
    li      t0, 1
    sll     t0, t0, s6
    srli    t1, t0, 1
    div     t2, s4, t0
    j       vbi_stage_loop

vbi_done:
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
# void fft_1d_vectorized(float *xr, float *xi, float *tr, float *ti, int n)
# Wrapper: bit_reverse then iterative butterfly (both vectorized)
# =============================================================================
    .globl fft_1d_vectorized
    .type  fft_1d_vectorized, @function
fft_1d_vectorized:
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

    # Step 1: vectorized bit reversal
    mv      a0, s0
    mv      a1, s1
    mv      a2, s4
    call    vec_bit_reverse_array

    # Step 2: vectorized iterative butterfly
    mv      a0, s0
    mv      a1, s1
    mv      a2, s2
    mv      a3, s3
    mv      a4, s4
    call    vec_butterfly_iterative

    ld      s4,   0(sp)
    ld      s3,   8(sp)
    ld      s2,  16(sp)
    ld      s1,  24(sp)
    ld      s0,  32(sp)
    ld      ra,  40(sp)
    addi    sp, sp, 48
    ret
