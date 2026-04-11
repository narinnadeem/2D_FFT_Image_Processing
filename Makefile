CC      = riscv64-linux-gnu-gcc
AS_FLAGS = -march=rv64gcv
CFLAGS  = -O0 -g -static -march=rv64gcv
QEMU    = qemu-riscv64
QFLAGS  = -cpu rv64,v=true,vlen=128

TARGET  = fft_m3

all: $(TARGET)

$(TARGET): main_m3.c fft_vectorized.s fft.s
	$(CC) $(CFLAGS) -o $(TARGET) main_m3.c fft_vectorized.s fft.s -lm

run: $(TARGET)
	$(QEMU) $(QFLAGS) ./$(TARGET)

clean:
	rm -f $(TARGET) *.o
