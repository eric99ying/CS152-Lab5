#include "mt-matmul.h"

void matmul_opt(unsigned int coreid, unsigned int ncores, const size_t lda,  const data_t* A, const data_t* B, data_t* C)
{
    int i, j, k;
    if (coreid % 2 == 0) {
        for (i = 0; i < 16; i += 2) {
            for (k = 0; k < 32; k += 4) {
                for (j = 0; j < 32; j += 2) {

                    data_t B1 = B[(k*lda) + j];
                    data_t B2 = B[(k*lda) + j + 1];
                    data_t B3 = B[((k + 1)*lda) + j];
                    data_t B4 = B[((k + 1)*lda) + j + 1];
                    data_t B5 = B[((k + 2)*lda) + j];
                    data_t B6 = B[((k + 2)*lda) + j + 1];
                    data_t B7 = B[((k + 3)*lda) + j];
                    data_t B8 = B[((k + 3)*lda) + j + 1];

                    data_t sum1 = 0; 
                    data_t sum2 = 0; 
                    data_t sum3 = 0; 
                    data_t sum4 = 0; 

                    if (k != 0) {
                        sum1 = C[i*lda + j];
                        sum2 = C[i*lda + j + 1];
                        sum3 = C[(i + 1)*lda + j];
                        sum4 = C[(i + 1)*lda + j + 1];
                    }

                    data_t A1 = A[(i*lda) + k];
                    data_t A2 = A[((i + 1)*lda) + k];
                    sum1 += A1 * B1;
                    sum2 += A1 * B2;
                    sum3 += A2 * B1;
                    sum4 += A2 * B2;

                    data_t A3 = A[(i*lda) + k + 1];
                    data_t A4 = A[((i + 1)*lda) + k + 1];
                    sum1 += A3 * B3;
                    sum2 += A3 * B4;
                    sum3 += A4 * B3;
                    sum4 += A4 * B4;    

                    data_t A5 = A[(i*lda) + k + 2];
                    data_t A6 = A[((i + 1)*lda) + k + 2];
                    sum1 += A5 * B5;
                    sum2 += A5 * B6;
                    sum3 += A6 * B5;
                    sum4 += A6 * B6;
                    
                    data_t A7 = A[(i*lda) + k + 3];
                    data_t A8 = A[((i + 1)*lda) + k + 3];
                    sum1 += A7 * B7;
                    sum2 += A7 * B8;
                    sum3 += A8 * B7;
                    sum4 += A8 * B8;

                    C[i*lda + j] = sum1;
                    C[i*lda + j + 1] = sum2;
                    C[(i + 1)*lda + j] = sum3;
                    C[(i + 1)*lda + j + 1] = sum4;
                }
            }
        }
    } else {
        for (i = 30; i >= 16; i -= 2) {
            for (k = 28; k >= 0; k -= 4) {
                for (j = 30; j >= 0; j -= 2) {

                    data_t B1 = B[(k*lda) + j];
                    data_t B2 = B[(k*lda) + j + 1];
                    data_t B3 = B[((k + 1)*lda) + j];
                    data_t B4 = B[((k + 1)*lda) + j + 1];
                    data_t B5 = B[((k + 2)*lda) + j];
                    data_t B6 = B[((k + 2)*lda) + j + 1];
                    data_t B7 = B[((k + 3)*lda) + j];
                    data_t B8 = B[((k + 3)*lda) + j + 1];

                    data_t sum1 = 0; 
                    data_t sum2 = 0; 
                    data_t sum3 = 0; 
                    data_t sum4 = 0; 

                    if (k != 28) {
                        sum1 = C[i*lda + j];
                        sum2 = C[i*lda + j + 1];
                        sum3 = C[(i + 1)*lda + j];
                        sum4 = C[(i + 1)*lda + j + 1];
                    }

                    data_t A1 = A[(i*lda) + k];
                    data_t A2 = A[((i + 1)*lda) + k];
                    sum1 += A1 * B1;
                    sum2 += A1 * B2;
                    sum3 += A2 * B1;
                    sum4 += A2 * B2;

                    data_t A3 = A[(i*lda) + k + 1];
                    data_t A4 = A[((i + 1)*lda) + k + 1];
                    sum1 += A3 * B3;
                    sum2 += A3 * B4;
                    sum3 += A4 * B3;
                    sum4 += A4 * B4;    

                    data_t A5 = A[(i*lda) + k + 2];
                    data_t A6 = A[((i + 1)*lda) + k + 2];
                    sum1 += A5 * B5;
                    sum2 += A5 * B6;
                    sum3 += A6 * B5;
                    sum4 += A6 * B6;
                    
                    data_t A7 = A[(i*lda) + k + 3];
                    data_t A8 = A[((i + 1)*lda) + k + 3];
                    sum1 += A7 * B7;
                    sum2 += A7 * B8;
                    sum3 += A8 * B7;
                    sum4 += A8 * B8;

                    C[i*lda + j] = sum1;
                    C[i*lda + j + 1] = sum2;
                    C[(i + 1)*lda + j] = sum3;
                    C[(i + 1)*lda + j + 1] = sum4;
                }
            }
        }
    }
}

