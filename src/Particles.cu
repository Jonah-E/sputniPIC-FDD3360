#include "Alloc.h"
#include "Particles.h"
#include "utils.h"
#include <cuda.h>
#include <cuda_runtime.h>

/** allocate particle arrays */
void particle_allocate(struct parameters* param, struct particles* part, int is)
{

    // set species ID
    part->species_ID = is;
    // number of particles
    part->nop = param->np[is];
    // maximum number of particles
    part->npmax = param->npMax[is];

    // choose a different number of mover iterations for ions and electrons
    if (param->qom[is] < 0) { // electrons
        part->NiterMover = param->NiterMover;
        part->n_sub_cycles = param->n_sub_cycles;
    } else { // ions: only one iteration
        part->NiterMover = 1;
        part->n_sub_cycles = 1;
    }

    // particles per cell
    part->npcelx = param->npcelx[is];
    part->npcely = param->npcely[is];
    part->npcelz = param->npcelz[is];
    part->npcel = part->npcelx * part->npcely * part->npcelz;

    // cast it to required precision
    part->qom = (FPpart) param->qom[is];

    long npmax = part->npmax;

    // initialize drift and thermal velocities
    // drift
    part->u0 = (FPpart) param->u0[is];
    part->v0 = (FPpart) param->v0[is];
    part->w0 = (FPpart) param->w0[is];
    // thermal
    part->uth = (FPpart) param->uth[is];
    part->vth = (FPpart) param->vth[is];
    part->wth = (FPpart) param->wth[is];

    //////////////////////////////
    /// ALLOCATION PARTICLE ARRAYS
    //////////////////////////////
    part->x = new FPpart[npmax];
    part->y = new FPpart[npmax];
    part->z = new FPpart[npmax];
    // allocate velocity
    part->u = new FPpart[npmax];
    part->v = new FPpart[npmax];
    part->w = new FPpart[npmax];
    // allocate charge = q * statistical weight
    part->q = new FPinterp[npmax];
}

/** deallocate */
void particle_deallocate(struct particles* part)
{
    // deallocate particle variables
    delete[] part->x;
    delete[] part->y;
    delete[] part->z;
    delete[] part->u;
    delete[] part->v;
    delete[] part->w;
    delete[] part->q;
}

/** particle mover */
int mover_PC(struct particles* part, struct EMfield* field, struct grid* grd,
             struct parameters* param)
{
    // print species and subcycling
    std::cout << "***  MOVER with SUBCYCLYING " << param->n_sub_cycles
              << " - species " << part->species_ID << " ***" << std::endl;

    // auxiliary variables
    FPpart dt_sub_cycling = (FPpart) param->dt / ((double) part->n_sub_cycles);
    FPpart dto2 = .5 * dt_sub_cycling, qomdt2 = part->qom * dto2 / param->c;
    FPpart omdtsq, denom, ut, vt, wt, udotb;

    // local (to the particle) electric and magnetic field
    FPfield Exl = 0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0, Bzl = 0.0;

    // interpolation densities
    int ix, iy, iz;
    FPfield weight[2][2][2];
    FPfield xi[2], eta[2], zeta[2];

    // intermediate particle position and velocity
    FPpart xptilde, yptilde, zptilde, uptilde, vptilde, wptilde;

    // start subcycling
    for (int i_sub = 0; i_sub < part->n_sub_cycles; i_sub++) {
        // move each particle with new fields
        for (int i = 0; i < part->nop; i++) {
            xptilde = part->x[i];
            yptilde = part->y[i];
            zptilde = part->z[i];
            // calculate the average velocity iteratively
            for (int innter = 0; innter < part->NiterMover; innter++) {
                // interpolation G-->P
                ix = 2 + int((part->x[i] - grd->xStart) * grd->invdx);
                iy = 2 + int((part->y[i] - grd->yStart) * grd->invdy);
                iz = 2 + int((part->z[i] - grd->zStart) * grd->invdz);

                // calculate weights
                xi[0] = part->x[i] - grd->XN[ix - 1][iy][iz];
                eta[0] = part->y[i] - grd->YN[ix][iy - 1][iz];
                zeta[0] = part->z[i] - grd->ZN[ix][iy][iz - 1];
                xi[1] = grd->XN[ix][iy][iz] - part->x[i];
                eta[1] = grd->YN[ix][iy][iz] - part->y[i];
                zeta[1] = grd->ZN[ix][iy][iz] - part->z[i];
                for (int ii = 0; ii < 2; ii++)
                    for (int jj = 0; jj < 2; jj++)
                        for (int kk = 0; kk < 2; kk++)
                            weight[ii][jj][kk] =
                                xi[ii] * eta[jj] * zeta[kk] * grd->invVOL;

                // set to zero local electric and magnetic field
                Exl = 0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0,
                Bzl = 0.0;

                for (int ii = 0; ii < 2; ii++)
                    for (int jj = 0; jj < 2; jj++)
                        for (int kk = 0; kk < 2; kk++) {
                            Exl += weight[ii][jj][kk] *
                                   field->Ex[ix - ii][iy - jj][iz - kk];
                            Eyl += weight[ii][jj][kk] *
                                   field->Ey[ix - ii][iy - jj][iz - kk];
                            Ezl += weight[ii][jj][kk] *
                                   field->Ez[ix - ii][iy - jj][iz - kk];
                            Bxl += weight[ii][jj][kk] *
                                   field->Bxn[ix - ii][iy - jj][iz - kk];
                            Byl += weight[ii][jj][kk] *
                                   field->Byn[ix - ii][iy - jj][iz - kk];
                            Bzl += weight[ii][jj][kk] *
                                   field->Bzn[ix - ii][iy - jj][iz - kk];
                        }

                // end interpolation
                omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
                denom = 1.0 / (1.0 + omdtsq);
                // solve the position equation
                ut = part->u[i] + qomdt2 * Exl;
                vt = part->v[i] + qomdt2 * Eyl;
                wt = part->w[i] + qomdt2 * Ezl;
                udotb = ut * Bxl + vt * Byl + wt * Bzl;
                // solve the velocity equation
                uptilde = (ut + qomdt2 * (vt * Bzl - wt * Byl +
                                          qomdt2 * udotb * Bxl)) *
                          denom;
                vptilde = (vt + qomdt2 * (wt * Bxl - ut * Bzl +
                                          qomdt2 * udotb * Byl)) *
                          denom;
                wptilde = (wt + qomdt2 * (ut * Byl - vt * Bxl +
                                          qomdt2 * udotb * Bzl)) *
                          denom;
                // update position
                part->x[i] = xptilde + uptilde * dto2;
                part->y[i] = yptilde + vptilde * dto2;
                part->z[i] = zptilde + wptilde * dto2;

            } // end of iteration
            // update the final position and velocity
            part->u[i] = 2.0 * uptilde - part->u[i];
            part->v[i] = 2.0 * vptilde - part->v[i];
            part->w[i] = 2.0 * wptilde - part->w[i];
            part->x[i] = xptilde + uptilde * dt_sub_cycling;
            part->y[i] = yptilde + vptilde * dt_sub_cycling;
            part->z[i] = zptilde + wptilde * dt_sub_cycling;

            //////////
            //////////
            ////////// BC

            // X-DIRECTION: BC particles
            if (part->x[i] > grd->Lx) {
                if (param->PERIODICX == true) { // PERIODIC
                    part->x[i] = part->x[i] - grd->Lx;
                } else { // REFLECTING BC
                    part->u[i] = -part->u[i];
                    part->x[i] = 2 * grd->Lx - part->x[i];
                }
            }

            if (part->x[i] < 0) {
                if (param->PERIODICX == true) { // PERIODIC
                    part->x[i] = part->x[i] + grd->Lx;
                } else { // REFLECTING BC
                    part->u[i] = -part->u[i];
                    part->x[i] = -part->x[i];
                }
            }

            // Y-DIRECTION: BC particles
            if (part->y[i] > grd->Ly) {
                if (param->PERIODICY == true) { // PERIODIC
                    part->y[i] = part->y[i] - grd->Ly;
                } else { // REFLECTING BC
                    part->v[i] = -part->v[i];
                    part->y[i] = 2 * grd->Ly - part->y[i];
                }
            }

            if (part->y[i] < 0) {
                if (param->PERIODICY == true) { // PERIODIC
                    part->y[i] = part->y[i] + grd->Ly;
                } else { // REFLECTING BC
                    part->v[i] = -part->v[i];
                    part->y[i] = -part->y[i];
                }
            }

            // Z-DIRECTION: BC particles
            if (part->z[i] > grd->Lz) {
                if (param->PERIODICZ == true) { // PERIODIC
                    part->z[i] = part->z[i] - grd->Lz;
                } else { // REFLECTING BC
                    part->w[i] = -part->w[i];
                    part->z[i] = 2 * grd->Lz - part->z[i];
                }
            }

            if (part->z[i] < 0) {
                if (param->PERIODICZ == true) { // PERIODIC
                    part->z[i] = part->z[i] + grd->Lz;
                } else { // REFLECTING BC
                    part->w[i] = -part->w[i];
                    part->z[i] = -part->z[i];
                }
            }

        } // end of subcycling
    }     // end of one particle

    return (0); // exit succcesfully
} // end of the mover

/** Interpolation Particle --> Grid: This is for species */
void interpP2G(struct particles* part, struct interpDensSpecies* ids,
               struct grid* grd)
{

    // arrays needed for interpolation
    FPpart weight[2][2][2];
    FPpart temp[2][2][2];
    FPpart xi[2], eta[2], zeta[2];

    // index of the cell
    int ix, iy, iz;

    for (register long long i = 0; i < part->nop; i++) {

        // determine cell: can we change to int()? is it faster?
        ix = 2 + int(floor((part->x[i] - grd->xStart) * grd->invdx));
        iy = 2 + int(floor((part->y[i] - grd->yStart) * grd->invdy));
        iz = 2 + int(floor((part->z[i] - grd->zStart) * grd->invdz));

        // distances from node
        xi[0] = part->x[i] - grd->XN[ix - 1][iy][iz];
        eta[0] = part->y[i] - grd->YN[ix][iy - 1][iz];
        zeta[0] = part->z[i] - grd->ZN[ix][iy][iz - 1];
        xi[1] = grd->XN[ix][iy][iz] - part->x[i];
        eta[1] = grd->YN[ix][iy][iz] - part->y[i];
        zeta[1] = grd->ZN[ix][iy][iz] - part->z[i];

        // calculate the weights for different nodes
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    weight[ii][jj][kk] =
                        part->q[i] * xi[ii] * eta[jj] * zeta[kk] * grd->invVOL;

        //////////////////////////
        // add charge density
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->rhon[ix - ii][iy - jj][iz - kk] +=
                        weight[ii][jj][kk] * grd->invVOL;

        ////////////////////////////
        // add current density - Jx
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] = part->u[i] * weight[ii][jj][kk];

        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->Jx[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        ////////////////////////////
        // add current density - Jy
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] = part->v[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->Jy[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        ////////////////////////////
        // add current density - Jz
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] = part->w[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->Jz[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        ////////////////////////////
        // add pressure pxx
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] =
                        part->u[i] * part->u[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->pxx[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        ////////////////////////////
        // add pressure pxy
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] =
                        part->u[i] * part->v[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->pxy[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        /////////////////////////////
        // add pressure pxz
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] =
                        part->u[i] * part->w[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->pxz[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        /////////////////////////////
        // add pressure pyy
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] =
                        part->v[i] * part->v[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->pyy[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        /////////////////////////////
        // add pressure pyz
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] =
                        part->v[i] * part->w[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->pyz[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;

        /////////////////////////////
        // add pressure pzz
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    temp[ii][jj][kk] =
                        part->w[i] * part->w[i] * weight[ii][jj][kk];
        for (int ii = 0; ii < 2; ii++)
            for (int jj = 0; jj < 2; jj++)
                for (int kk = 0; kk < 2; kk++)
                    ids->pzz[ix - ii][iy - jj][iz - kk] +=
                        temp[ii][jj][kk] * grd->invVOL;
    }
}

// x*(sizeY + sizeZ) + y*sizeZ + z
__global__ void move_particle_kernel(struct particles_gpu part,
                                     struct grid_gpu grd,
                                     struct EMfield_gpu field,
                                     struct parameters_gpu param)
{
    const unsigned int tid = threadIdx.x + blockDim.x * blockIdx.x;

    if (tid >= part.nop){
        return;
    }

    long idx, idy;

    // auxiliary variables
    FPpart dt_sub_cycling = (FPpart) param.dt / ((double) part.n_sub_cycles);
    FPpart dto2 = .5 * dt_sub_cycling, qomdt2 = part.qom * dto2 / param.c;
    FPpart omdtsq, denom, ut, vt, wt, udotb;

    // local (to the particle) electric and magnetic field
    FPfield Exl = 0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0, Bzl = 0.0;

    // interpolation densities
    int ix, iy, iz;
    FPfield weight[2 * 2 * 2];
    FPfield xi[2], eta[2], zeta[2];

    // intermediate particle position and velocity
    FPpart xptilde, yptilde, zptilde, uptilde, vptilde, wptilde;

    for (int i_sub = 0; i_sub < part.n_sub_cycles; i_sub++) {
        xptilde = part.x[tid];
        yptilde = part.y[tid];
        zptilde = part.z[tid];
        // calculate the average velocity iteratively
        for (int innter = 0; innter < part.NiterMover; innter++) {
            // interpolation G-.P
            ix = 2 + int((part.x[tid] - grd.xStart) * grd.invdx);
            iy = 2 + int((part.y[tid] - grd.yStart) * grd.invdy);
            iz = 2 + int((part.z[tid] - grd.zStart) * grd.invdz);

            // calculate weights
            idx = get_idx(ix - 1, iy, iz, grd.nyn, grd.nzn);
            xi[0] = part.x[tid] - grd.XN_flat[idx];
            idx = get_idx(ix, iy - 1, iz, grd.nyn, grd.nzn);
            eta[0] = part.y[tid] - grd.YN_flat[idx];
            idx = get_idx(ix, iy, iz - 1, grd.nyn, grd.nzn);
            zeta[0] = part.z[tid] - grd.ZN_flat[idx];

            idx = get_idx(ix, iy, iz, grd.nyn, grd.nzn);
            xi[1] = grd.XN_flat[idx] - part.x[tid];
            eta[1] = grd.YN_flat[idx] - part.y[tid];
            zeta[1] = grd.ZN_flat[idx] - part.z[tid];

            for (int ii = 0; ii < 2; ii++)
                for (int jj = 0; jj < 2; jj++)
                    for (int kk = 0; kk < 2; kk++) {
                        idx = get_idx(ii, jj, kk, 2, 2);
                        weight[idx] = xi[ii] * eta[jj] * zeta[kk] * grd.invVOL;
                    }

            // set to zero local electric and magnetic field
            Exl = 0.0, Eyl = 0.0, Ezl = 0.0, Bxl = 0.0, Byl = 0.0, Bzl = 0.0;

            for (int ii = 0; ii < 2; ii++)
                for (int jj = 0; jj < 2; jj++)
                    for (int kk = 0; kk < 2; kk++) {
                        idx = get_idx(ii, jj, kk, 2, 2);
                        idy = get_idx(ix - ii, iy - jj, iz - kk, grd.nyn,
                                      grd.nzn);
                        Exl += weight[idx] * field.Ex_flat[idy];
                        Eyl += weight[idx] * field.Ey_flat[idy];
                        Ezl += weight[idx] * field.Ez_flat[idy];
                        Bxl += weight[idx] * field.Bxn_flat[idy];
                        Byl += weight[idx] * field.Byn_flat[idy];
                        Bzl += weight[idx] * field.Bzn_flat[idy];
                    }

            // end interpolation
            omdtsq = qomdt2 * qomdt2 * (Bxl * Bxl + Byl * Byl + Bzl * Bzl);
            denom = 1.0 / (1.0 + omdtsq);
            // solve the position equation
            ut = part.u[tid] + qomdt2 * Exl;
            vt = part.v[tid] + qomdt2 * Eyl;
            wt = part.w[tid] + qomdt2 * Ezl;
            udotb = ut * Bxl + vt * Byl + wt * Bzl;
            // solve the velocity equation
            uptilde =
                (ut + qomdt2 * (vt * Bzl - wt * Byl + qomdt2 * udotb * Bxl)) *
                denom;
            vptilde =
                (vt + qomdt2 * (wt * Bxl - ut * Bzl + qomdt2 * udotb * Byl)) *
                denom;
            wptilde =
                (wt + qomdt2 * (ut * Byl - vt * Bxl + qomdt2 * udotb * Bzl)) *
                denom;
            // update position
            part.x[tid] = xptilde + uptilde * dto2;
            part.y[tid] = yptilde + vptilde * dto2;
            part.z[tid] = zptilde + wptilde * dto2;

        } // end of iteration
        // update the final position and velocity
        part.u[tid] = 2.0 * uptilde - part.u[tid];
        part.v[tid] = 2.0 * vptilde - part.v[tid];
        part.w[tid] = 2.0 * wptilde - part.w[tid];
        part.x[tid] = xptilde + uptilde * dt_sub_cycling;
        part.y[tid] = yptilde + vptilde * dt_sub_cycling;
        part.z[tid] = zptilde + wptilde * dt_sub_cycling;

        //////////
        //////////
        ////////// BC

        // X-DIRECTION: BC particles
        if (part.x[tid] > grd.Lx) {
            if (param.PERIODICX == true) { // PERIODIC
                part.x[tid] = part.x[tid] - grd.Lx;
            } else { // REFLECTING BC
                part.u[tid] = -part.u[tid];
                part.x[tid] = 2 * grd.Lx - part.x[tid];
            }
        }

        if (part.x[tid] < 0) {
            if (param.PERIODICX == true) { // PERIODIC
                part.x[tid] = part.x[tid] + grd.Lx;
            } else { // REFLECTING BC
                part.u[tid] = -part.u[tid];
                part.x[tid] = -part.x[tid];
            }
        }

        // Y-DIRECTION: BC particles
        if (part.y[tid] > grd.Ly) {
            if (param.PERIODICY == true) { // PERIODIC
                part.y[tid] = part.y[tid] - grd.Ly;
            } else { // REFLECTING BC
                part.v[tid] = -part.v[tid];
                part.y[tid] = 2 * grd.Ly - part.y[tid];
            }
        }

        if (part.y[tid] < 0) {
            if (param.PERIODICY == true) { // PERIODIC
                part.y[tid] = part.y[tid] + grd.Ly;
            } else { // REFLECTING BC
                part.v[tid] = -part.v[tid];
                part.y[tid] = -part.y[tid];
            }
        }

        // Z-DIRECTION: BC particles
        if (part.z[tid] > grd.Lz) {
            if (param.PERIODICZ == true) { // PERIODIC
                part.z[tid] = part.z[tid] - grd.Lz;
            } else { // REFLECTING BC
                part.w[tid] = -part.w[tid];
                part.z[tid] = 2 * grd.Lz - part.z[tid];
            }
        }

        if (part.z[tid] < 0) {
            if (param.PERIODICZ == true) { // PERIODIC
                part.z[tid] = part.z[tid] + grd.Lz;
            } else { // REFLECTING BC
                part.w[tid] = -part.w[tid];
                part.z[tid] = -part.z[tid];
            }
        }

    } // end of subcycling
}

void particles_cpy(struct particles_gpu* dst, struct particles* src)
{
    dst->NiterMover = src->NiterMover;
    dst->n_sub_cycles = src->n_sub_cycles;
    dst->qom = src->qom;

    if (src->nop <= src->npmax) {
        dst->nop = src->nop;
    } else {
        dst->nop = src->npmax;
        printf("Error: Number of particles is larger than max (%ld > %ld), using "
               "max for copy",
               src->nop, src->npmax);
    }
}
void particles_cpy_to_gpu(struct particles_gpu* dst, struct particles* src)
{

    cudaError_t deviceError;
    deviceError = cudaMemcpy(dst->x, src->x, sizeof(FPpart) * dst->nop,
                             cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->y, src->y, sizeof(FPpart) * dst->nop,
                             cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->z, src->z, sizeof(FPpart) * dst->nop,
                             cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }

    deviceError = cudaMemcpy(dst->u, src->u, sizeof(FPpart) * dst->nop,
                             cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->v, src->v, sizeof(FPpart) * dst->nop,
                             cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->w, src->w, sizeof(FPpart) * dst->nop,
                             cudaMemcpyHostToDevice);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
}

void particles_cpy_from_gpu(struct particles* dst, struct particles_gpu* src)
{

    cudaError_t deviceError;
    deviceError = cudaMemcpy(dst->x, src->x, sizeof(FPpart) * src->nop,
                             cudaMemcpyDeviceToHost);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->y, src->y, sizeof(FPpart) * src->nop,
                             cudaMemcpyDeviceToHost);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->z, src->z, sizeof(FPpart) * src->nop,
                             cudaMemcpyDeviceToHost);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }

    deviceError = cudaMemcpy(dst->u, src->u, sizeof(FPpart) * src->nop,
                             cudaMemcpyDeviceToHost);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->v, src->v, sizeof(FPpart) * src->nop,
                             cudaMemcpyDeviceToHost);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
    deviceError = cudaMemcpy(dst->w, src->w, sizeof(FPpart) * src->nop,
                             cudaMemcpyDeviceToHost);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
    }
}

cudaError_t particle_allocate_gpu(struct particles_gpu* gpu_part, long npmax)
{
    cudaError_t deviceError;
    deviceError = cudaMalloc(&gpu_part->x, sizeof(FPpart) * npmax);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError = cudaMalloc(&gpu_part->y, sizeof(FPpart) * npmax);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError = cudaMalloc(&gpu_part->z, sizeof(FPpart) * npmax);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }

    deviceError = cudaMalloc(&gpu_part->u, sizeof(FPpart) * npmax);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError = cudaMalloc(&gpu_part->v, sizeof(FPpart) * npmax);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }
    deviceError = cudaMalloc(&gpu_part->w, sizeof(FPpart) * npmax);
    if (deviceError != cudaSuccess) {
        printCudaError(deviceError);
        return deviceError;
    }

    return cudaSuccess;
}

void particle_deallocate_gpu(struct particles_gpu* gpu_part)
{
    cudaFree(gpu_part->x);
    cudaFree(gpu_part->y);
    cudaFree(gpu_part->z);

    cudaFree(gpu_part->u);
    cudaFree(gpu_part->v);
    cudaFree(gpu_part->w);
}

void deallocate_gpu(struct particles_gpu* gpu_part, struct grid_gpu* gpu_grd,
                    struct EMfield_gpu* gpu_field)
{
    particle_deallocate_gpu(gpu_part);

    grid_deallocate_gpu(gpu_grd);

    emfield_deallocate_gpu(gpu_field);
}

#define MAX_THREADS_PER_BLOCK 1024

/** particle mover */
int mover_PC_gpu(struct particles* part, struct EMfield* field,
                 struct grid* grd, struct parameters* param)
{
    // print species and subcycling
    std::cout << "***  MOVER with SUBCYCLYING " << param->n_sub_cycles
              << " - species " << part->species_ID << " ***" << std::endl;

    struct parameters_gpu gpu_param;
    struct particles_gpu gpu_part;
    struct grid_gpu gpu_grd;
    struct EMfield_gpu gpu_field;

    long grd_arrays_size = grd->nxn * grd->nyn * grd->nzn;
    cudaError_t deviceError;

    deviceError = particle_allocate_gpu(&gpu_part, part->npmax);
    if (deviceError != cudaSuccess) {
        particle_deallocate_gpu(&gpu_part);
        return 1;
    }

    deviceError = grid_allocate_gpu(&gpu_grd, grd_arrays_size);
    if (deviceError != cudaSuccess) {
        grid_deallocate_gpu(&gpu_grd);
        particle_deallocate_gpu(&gpu_part);
        return 1;
    }

    deviceError = emfield_allocate_gpu(&gpu_field, grd_arrays_size);
    if (deviceError != cudaSuccess) {
        emfield_deallocate_gpu(&gpu_field);
        grid_deallocate_gpu(&gpu_grd);
        particle_deallocate_gpu(&gpu_part);
        return 1;
    }

    particles_cpy(&gpu_part, part);
    parameters_cpy(&gpu_param, param);
    grid_cpy(&gpu_grd, grd);

    particles_cpy_to_gpu(&gpu_part, part);
    grid_cpy_to_gpu(&gpu_grd, grd);
    emfield_cpy_to_gpu(&gpu_field, field, grd_arrays_size);

    dim3 grid((gpu_part.nop + MAX_THREADS_PER_BLOCK - 1) /
                  MAX_THREADS_PER_BLOCK,
              1, 1);
    dim3 block(MAX_THREADS_PER_BLOCK, 1, 1);

    move_particle_kernel<<<grid, block>>>(gpu_part, gpu_grd, gpu_field,
                                          gpu_param);
    deviceError = cudaGetLastError();
    if (deviceError != cudaSuccess) {
        printf("Error when running GPU: %s (%d)\n",
               cudaGetErrorString(deviceError), deviceError);
    }

    particles_cpy_from_gpu(part, &gpu_part);

    deallocate_gpu(&gpu_part, &gpu_grd, &gpu_field);
    return (0); // exit succcesfully
} // end of one particle
