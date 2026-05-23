CC      = riscv64-linux-gnu-gcc
CFLAGS  = -O0 -g -static -march=rv64gcv
QEMU    = qemu-riscv64
QFLAGS  = -cpu rv64,v=true,vlen=128
HCC     = cc

TARGET  = fft_m4
VIZ     = visualize

all: $(TARGET) $(VIZ)

$(TARGET): main_m4.c fft_vectorized.s fft.s
	$(CC) $(CFLAGS) -o $(TARGET) main_m4.c fft_vectorized.s fft.s -lm

$(VIZ): visualize.c
	$(HCC) visualize.c -o $(VIZ) -lm

run: $(TARGET)
	$(QEMU) $(QFLAGS) ./$(TARGET)

viz: $(VIZ)
	./$(VIZ)

clean:
	rm -f $(TARGET) $(VIZ) fft_result.bin fft_output.png edges_output.png
