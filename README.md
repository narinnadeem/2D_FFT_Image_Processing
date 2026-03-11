# 2D FFT Image Processing (Milestone 1)

This project implements a 2D Fast Fourier Transform (FFT) pipeline in C for image edge detection. It is designed to be portable for future RISC-V assembly implementation.

## Project Structure
- `main.c`: Orchestrates image loading, FFT processing, and saving results.
- `fft_implementation.c/h`: Core 1D and 2D FFT logic, including bit-reversal and butterfly computations.
- `math.c/h`: Custom Taylor Series approximations for `sin`, `cos`, and `exp`.
- `stb_image.h` / `stb_image_write.h`: External libraries for image handling.

## How to Compile and Run
Ensure you have a C compiler (like `gcc`) and the test images (`8x8.png`, `jet.png`) in the same folder.

### Run the Image Processor:
gcc main.c fft_implementation.c math.c -o fft_proc -lm
./fft_proc

### Run the Test Suite:
gcc FFT-milestone1-testsuite.c fft_implementation.c math.c -o test_suite -lm
./test_suite