#include "mt-vvadd.h"

void vvadd_opt(int coreid, int ncores, size_t n, const data_t* x, const data_t* y, data_t* z)
{
    size_t i;
    // Interleave accesses
    int block = n / ncores;
    int stop;
    if ((coreid + 1) * block > n) {
        stop = n;
    } else {
        stop = (coreid + 1) * block;
    }

    for (i = coreid * block; i < stop; i += 1) {
        z[i] = x[i] + y[i];
    }
}
