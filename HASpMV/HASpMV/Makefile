
CC = gcc
CFLAGS = -O3 -I ../include -fopenmp -march=native -DCORE_P=16 -DCORE_E=8

SOURCES = haspmv.c naive.c partition.c reorder.c avx2.c loop.c
TARGETS = $(SOURCES:.c=)

all: $(TARGETS)

%: %.c
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(TARGETS)
