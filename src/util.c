#include <stdlib.h>
#include <time.h>
#include "util.h"

#ifdef _WIN32
#include <windows.h>
#endif

// --- random utilities ---
double unif01(void){ return (double)rand() / (double)RAND_MAX; }
int max_int(int a, int b){ return (a>b)?a:b; }
int min_int(int a, int b){ return (a<b)?a:b; }

// --- high-resolution timer (for measuring runtime) ---
double now_sec(void){
#ifdef _WIN32
    LARGE_INTEGER freq, t;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t);
    return (double)t.QuadPart / (double)freq.QuadPart;
#elif defined(CLOCK_MONOTONIC)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
#else
    return (double)clock() / (double)CLOCKS_PER_SEC;
#endif
}
