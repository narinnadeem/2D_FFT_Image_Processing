# 2D FFT Image Processing

## Milestone 1 — C Implementation (branch: main)
2D FFT pipeline in C for image edge detection using custom math functions.

### Compile & Run
gcc main.c fft_implementation.c math.c -o fft_proc -lm
./fft_proc

## Milestone 2 — RISC-V Assembly (branch: milestone-2)
1D FFT implemented in RISC-V assembly, emulated on QEMU.

### Files
- fft.s — RISC-V assembly: sin, cos, twiddle factors, bit reversal, iterative & recursive FFT
- main_m2.c — C driver for testing
- Makefile — Build and run

### Compile & Run
make
make run
