#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <cuda.h>

#define printCudaError(cuda_returned_error_code)                               \
    {                                                                          \
        cudaErrorPrint((cuda_returned_error_code), __FILE__, __LINE__);        \
    }

inline void cudaErrorPrint(cudaError_t code, const char* file, int line)
{
    fprintf(stderr, "CUDA Error: %s (%d) %s %d\n", cudaGetErrorString(code),
            code, file, line);
}

#endif
