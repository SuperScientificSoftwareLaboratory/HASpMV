# generate_makefile.py
import sys

input_value = sys.argv[1]

core_p = 16
core_e = 8

if input_value == '12':
    core_p = 16
    core_e = 8
elif input_value == '13':
    core_p = 16
    core_e = 16

makefile_content = f"""
CC = gcc
CFLAGS = -O3 -I ../include -fopenmp -march=native -DCORE_P={core_p} -DCORE_E={core_e}

SOURCES = haspmv.c naive.c partition.c reorder.c avx2.c loop.c
TARGETS = $(SOURCES:.c=)

all: $(TARGETS)

%: %.c
\t$(CC) $(CFLAGS) $< -o $@

clean:
\trm -f $(TARGETS)
"""

with open('Makefile', 'w') as makefile:
    makefile.write(makefile_content)
