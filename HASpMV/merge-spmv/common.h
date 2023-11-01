#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>

#include <mm_malloc.h>
#include <x86intrin.h>
#include <immintrin.h>
#include <nmmintrin.h>

#include <omp.h>

#ifndef MAT_VAL_TYPE
#define MAT_VAL_TYPE double
#endif

#ifndef BENCH_REPEAT
#define BENCH_REPEAT 1000
#endif

#ifndef MAT_PTR_TYPE
#define MAT_PTR_TYPE int
#endif