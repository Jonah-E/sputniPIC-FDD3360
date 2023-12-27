#ifndef EMFIELD_H
#define EMFIELD_H

#include "Alloc.h"
#include "Grid.h"

/** structure with field information */
struct EMfield {
    // field arrays: 4D arrays

    /* Electric field defined on nodes: last index is component */
    FPfield*** Ex;
    FPfield* Ex_flat;
    FPfield*** Ey;
    FPfield* Ey_flat;
    FPfield*** Ez;
    FPfield* Ez_flat;
    /* Magnetic field defined on nodes: last index is component */
    FPfield*** Bxn;
    FPfield* Bxn_flat;
    FPfield*** Byn;
    FPfield* Byn_flat;
    FPfield*** Bzn;
    FPfield* Bzn_flat;
};

/** allocate electric and magnetic field */
void field_allocate(struct grid*, struct EMfield*);

/** deallocate electric and magnetic field */
void field_deallocate(struct grid*, struct EMfield*);

/** GPU */
struct EMfield_gpu {
    FPfield* Ex_flat;
    FPfield* Ey_flat;
    FPfield* Ez_flat;
    FPfield* Bxn_flat;
    FPfield* Byn_flat;
    FPfield* Bzn_flat;
};

cudaError_t emfield_allocate_gpu(struct EMfield_gpu* gpu_field,
                                 size_t grd_arrays_size);
void emfield_deallocate_gpu(struct EMfield_gpu* gpu_field);
void emfield_cpy_to_gpu(struct EMfield_gpu* dst, struct EMfield* src,
                      size_t grd_arrays_size);

#endif
