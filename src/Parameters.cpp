#include "Parameters.h"

void parameters_cpy(struct parameters_gpu* dst, struct parameters* src)
{
    dst->dt = src->dt;
    dst->c = src->c;

    dst->PERIODICX = src->PERIODICX;
    dst->PERIODICY = src->PERIODICY;
    dst->PERIODICZ = src->PERIODICZ;
}

