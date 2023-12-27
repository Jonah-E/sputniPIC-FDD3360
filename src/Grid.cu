#include "Grid.h"
#include "utils.h"

cudaError_t grid_allocate_gpu(struct grid_gpu* gpu_grd, size_t grd_arrays_size)
{
    cudaError_t deviceError;
    deviceError =
        cudaMalloc(&gpu_grd->XN_flat, sizeof(FPfield) * grd_arrays_size);
    if(deviceError != cudaSuccess){
        printCudaError(deviceError);
        return deviceError;
    }

    deviceError =
        cudaMalloc(&gpu_grd->YN_flat, sizeof(FPfield) * grd_arrays_size);
    if(deviceError != cudaSuccess){
        printCudaError(deviceError);
        return deviceError;
    }

    deviceError =
        cudaMalloc(&gpu_grd->ZN_flat, sizeof(FPfield) * grd_arrays_size);
    if(deviceError != cudaSuccess){
        printCudaError(deviceError);
        return deviceError;
    }

    return cudaSuccess;
}

void grid_deallocate_gpu(struct grid_gpu* gpu_grd) {
    cudaFree(gpu_grd->XN_flat);
    cudaFree(gpu_grd->YN_flat);
    cudaFree(gpu_grd->ZN_flat);
}

void grid_cpy(struct grid_gpu* dst, struct grid* src)
{
    dst->nyn = src->nyn;
    dst->nzn = src->nzn;
    dst->invdx = src->invdx;
    dst->invdy = src->invdy;
    dst->invdz = src->invdz;
    dst->xStart = src->xStart;
    dst->yStart = src->yStart;
    dst->zStart = src->zStart;

    dst->invVOL = src->invVOL;

    dst->Lx = src->Lx;
    dst->Ly = src->Ly;
    dst->Lz = src->Lz;
}

void grid_cpy_to_gpu(struct grid_gpu* dst, struct grid* src)
{
    cudaError_t deviceError;
    long grd_arrays_size = src->nxn * src->nyn * src->nzn;

    deviceError = cudaMemcpy(dst->XN_flat, src->XN_flat,
                              sizeof(FPfield) * grd_arrays_size,
                              cudaMemcpyHostToDevice);
    if(deviceError != cudaSuccess){
        printCudaError(deviceError);
    }

    deviceError = cudaMemcpy(dst->YN_flat, src->YN_flat,
                              sizeof(FPfield) * grd_arrays_size,
                              cudaMemcpyHostToDevice);
    if(deviceError != cudaSuccess){
        printCudaError(deviceError);
    }

    deviceError = cudaMemcpy(dst->ZN_flat, src->ZN_flat,
                              sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if(deviceError != cudaSuccess){
        printCudaError(deviceError);
    }
}
