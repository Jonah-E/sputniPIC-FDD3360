#include "EMfield.h"
#include "utils.h"
#include <cuda.h>

cudaError_t emfield_allocate_gpu(struct EMfield_gpu* gpu_field,
                                 size_t grd_arrays_size)
{
    cudaError_t deviceError;
    deviceError =
        cudaMalloc(&gpu_field->Ex_flat, sizeof(FPfield) * grd_arrays_size);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError =
        cudaMalloc(&gpu_field->Ey_flat, sizeof(FPfield) * grd_arrays_size);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError =
        cudaMalloc(&gpu_field->Ez_flat, sizeof(FPfield) * grd_arrays_size);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError =
        cudaMalloc(&gpu_field->Bxn_flat, sizeof(FPfield) * grd_arrays_size);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError =
        cudaMalloc(&gpu_field->Byn_flat, sizeof(FPfield) * grd_arrays_size);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError =
        cudaMalloc(&gpu_field->Bzn_flat, sizeof(FPfield) * grd_arrays_size);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }

    return cudaSuccess;
}

void emfield_deallocate_gpu(struct EMfield_gpu* gpu_field)
{
    cudaFree(gpu_field->Ex_flat);
    cudaFree(gpu_field->Ey_flat);
    cudaFree(gpu_field->Ez_flat);

    cudaFree(gpu_field->Bxn_flat);
    cudaFree(gpu_field->Byn_flat);
    cudaFree(gpu_field->Bzn_flat);
}

void emfield_cpy_to_gpu(struct EMfield_gpu* dst, struct EMfield* src,
                      size_t grd_arrays_size)
{
    cudaError_t deviceError;
    deviceError =
        cudaMemcpy(dst->Ex_flat, src->Ex_flat,
                   sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError =
        cudaMemcpy(dst->Ey_flat, src->Ey_flat,
                   sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError =
        cudaMemcpy(dst->Ez_flat, src->Ez_flat,
                   sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }

    deviceError =
        cudaMemcpy(dst->Bxn_flat, src->Bxn_flat,
                   sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError =
        cudaMemcpy(dst->Byn_flat, src->Byn_flat,
                   sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError =
        cudaMemcpy(dst->Bzn_flat, src->Bzn_flat,
                   sizeof(FPfield) * grd_arrays_size, cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
}
