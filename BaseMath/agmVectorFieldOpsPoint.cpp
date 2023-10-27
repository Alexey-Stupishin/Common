#include "stdDefinitions.h"
#include <math.h>

#include "agmVectorFieldOps.h"
#include "agmScalarFieldOps.h"

#define fidx(kx, ky, kz) (ky)+(kz)*N[1]][(kx)

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::cross_plane(CagmVectorFieldOps *a, const CagmVectorFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    REALTYPE_A tx, ty, tz;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            double ax = a->fieldX[fidx(kx, ky, kz)];
            double ay = a->fieldY[fidx(kx, ky, kz)];
            double az = a->fieldZ[fidx(kx, ky, kz)];
            double bx = b->fieldX[fidx(kx, ky, kz)];
            double by = b->fieldY[fidx(kx, ky, kz)];
            double bz = b->fieldZ[fidx(kx, ky, kz)];

            tx = a->fieldY[fidx(kx, ky, kz)] * b->fieldZ[fidx(kx, ky, kz)] - a->fieldZ[fidx(kx, ky, kz)] * b->fieldY[fidx(kx, ky, kz)];
            ty = a->fieldZ[fidx(kx, ky, kz)] * b->fieldX[fidx(kx, ky, kz)] - a->fieldX[fidx(kx, ky, kz)] * b->fieldZ[fidx(kx, ky, kz)];
            tz = a->fieldX[fidx(kx, ky, kz)] * b->fieldY[fidx(kx, ky, kz)] - a->fieldY[fidx(kx, ky, kz)] * b->fieldX[fidx(kx, ky, kz)];
            fieldX[fidx(kx, ky, kz)] = tx;
            fieldY[fidx(kx, ky, kz)] = ty;
            fieldZ[fidx(kx, ky, kz)] = tz;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::rot_plane(CagmVectorFieldOps *a, int kz)
{
    REALTYPE_A zy, yz, xz, zx, yx, xy;
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
        {
            if (kx == 0)
            {
                zx = -3 * a->fieldZ[fidx(0, ky, kz)] + 4 * a->fieldZ[fidx(1, ky, kz)] - a->fieldZ[fidx(2, ky, kz)];
                yx = -3 * a->fieldY[fidx(0, ky, kz)] + 4 * a->fieldY[fidx(1, ky, kz)] - a->fieldY[fidx(2, ky, kz)];
            }
            else if (kx == N[0] - 1)
            {
                zx = a->fieldZ[fidx(N[0] - 3, ky, kz)] - 4 * a->fieldZ[fidx(N[0] - 2, ky, kz)] + 3 * a->fieldZ[fidx(N[0] - 1, ky, kz)];
                yx = a->fieldY[fidx(N[0] - 3, ky, kz)] - 4 * a->fieldY[fidx(N[0] - 2, ky, kz)] + 3 * a->fieldY[fidx(N[0] - 1, ky, kz)];
            }
            else
            {
                zx = a->fieldZ[fidx(kx + 1, ky, kz)] - a->fieldZ[fidx(kx - 1, ky, kz)];
                yx = a->fieldY[fidx(kx + 1, ky, kz)] - a->fieldY[fidx(kx - 1, ky, kz)];
            }

            if (ky == 0)
            {
                zy = -3 * a->fieldZ[fidx(kx, 0, kz)] + 4 * a->fieldZ[fidx(kx, 1, kz)] - a->fieldZ[fidx(kx, 2, kz)];
                xy = -3 * a->fieldX[fidx(kx, 0, kz)] + 4 * a->fieldX[fidx(kx, 1, kz)] - a->fieldX[fidx(kx, 2, kz)];
            }
            else if (ky == N[1] - 1)
            {
                zy = a->fieldZ[fidx(kx, N[1] - 3, kz)] - 4 * a->fieldZ[fidx(kx, N[1] - 2, kz)] + 3 * a->fieldZ[fidx(kx, N[1] - 1, kz)];
                xy = a->fieldX[fidx(kx, N[1] - 3, kz)] - 4 * a->fieldX[fidx(kx, N[1] - 2, kz)] + 3 * a->fieldX[fidx(kx, N[1] - 1, kz)];
            }
            else
            {
                zy = a->fieldZ[fidx(kx, ky + 1, kz)] - a->fieldZ[fidx(kx, ky - 1, kz)];
                xy = a->fieldX[fidx(kx, ky + 1, kz)] - a->fieldX[fidx(kx, ky - 1, kz)];
            }

            if (kz == 0)
            {
                xz = -3 * a->fieldX[fidx(kx, ky, 0)] + 4 * a->fieldX[fidx(kx, ky, 1)] - a->fieldX[fidx(kx, ky, 2)];
                yz = -3 * a->fieldY[fidx(kx, ky, 0)] + 4 * a->fieldY[fidx(kx, ky, 1)] - a->fieldY[fidx(kx, ky, 2)];
            }
            else if (kz == N[2] - 1)
            {
                xz = a->fieldX[fidx(kx, ky, N[2] - 3)] - 4 * a->fieldX[fidx(kx, ky, N[2] - 2)] + 3 * a->fieldX[fidx(kx, ky, N[2] - 1)];
                yz = a->fieldY[fidx(kx, ky, N[2] - 3)] - 4 * a->fieldY[fidx(kx, ky, N[2] - 2)] + 3 * a->fieldY[fidx(kx, ky, N[2] - 1)];
            }
            else
            {
                xz = a->fieldX[fidx(kx, ky, kz + 1)] - a->fieldX[fidx(kx, ky, kz - 1)];
                yz = a->fieldY[fidx(kx, ky, kz + 1)] - a->fieldY[fidx(kx, ky, kz - 1)];
            }

            fieldX[fidx(kx, ky, kz)] = (zy*step[1] - yz*step[2])*0.5;
            fieldY[fidx(kx, ky, kz)] = (xz*step[2] - zx*step[0])*0.5;
            fieldZ[fidx(kx, ky, kz)] = (yx*step[0] - xy*step[1])*0.5;
            //fieldX[fidx(kx, ky, kz)] = (zy - yz)*0.5;
            //fieldY[fidx(kx, ky, kz)] = (xz - zx)*0.5;
            //fieldZ[fidx(kx, ky, kz)] = (yx - xy)*0.5;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::grad_plane(CagmScalarFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    REALTYPE_A dx, dy, dz;
    for (int ky = 0; ky < N[1]; ky++)
        for (int kx = 0; kx < N[0]; kx++)
        {
            if (kx == 0)
                dx = -3 * a->field[fidx(0, ky, kz)] + 4 * a->field[fidx(1, ky, kz)] - a->field[fidx(2, ky, kz)];
            else if (kx == N[0] - 1)
                dx = a->field[fidx(N[0] - 3, ky, kz)] - 4 * a->field[fidx(N[0] - 2, ky, kz)] + 3 * a->field[fidx(N[0] - 1, ky, kz)];
            else
                dx = a->field[fidx(kx + 1, ky, kz)] - a->field[fidx(kx - 1, ky, kz)];

            if (ky == 0)
                dy = -3 * a->field[fidx(kx, 0, kz)] + 4 * a->field[fidx(kx, 1, kz)] - a->field[fidx(kx, 2, kz)];
            else if (ky == N[1] - 1)
                dy = a->field[fidx(kx, N[1] - 3, kz)] - 4 * a->field[fidx(kx, N[1] - 2, kz)] + 3 * a->field[fidx(kx, N[1] - 1, kz)];
            else
                dy = a->field[fidx(kx, ky + 1, kz)] - a->field[fidx(kx, ky - 1, kz)];

            if (kz == 0)
                dz = -3 * a->field[fidx(kx, ky, 0)] + 4 * a->field[fidx(kx, ky, 1)] - a->field[fidx(kx, ky, 2)];
            else if (kz == N[2] - 1)
                dz = a->field[fidx(kx, ky, N[2] - 3)] - 4 * a->field[fidx(kx, ky, N[2] - 2)] + 3 * a->field[fidx(kx, ky, N[2] - 1)];
            else
                dz = a->field[fidx(kx, ky, kz + 1)] - a->field[fidx(kx, ky, kz - 1)];

            fieldX[fidx(kx, ky, kz)] = dx*0.5*step[0];
            fieldY[fidx(kx, ky, kz)] = dy*0.5*step[1];
            fieldZ[fidx(kx, ky, kz)] = dz*0.5*step[2];
            //fieldX[fidx(kx, ky, kz)] = dx*0.5;
            //fieldY[fidx(kx, ky, kz)] = dy*0.5;
            //fieldZ[fidx(kx, ky, kz)] = dz*0.5;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult_plane(REALTYPE_A c, CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] * c;
            fieldY[fidx(kx, ky, kz)] = a->fieldY[fidx(kx, ky, kz)] * c;
            fieldZ[fidx(kx, ky, kz)] = a->fieldZ[fidx(kx, ky, kz)] * c;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::mult_plane(CagmScalarFieldOps *c, CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = (a->fieldX[fidx(kx, ky, kz)]) * (c->field[fidx(kx, ky, kz)]);
            fieldY[fidx(kx, ky, kz)] = (a->fieldY[fidx(kx, ky, kz)]) * (c->field[fidx(kx, ky, kz)]);
            fieldZ[fidx(kx, ky, kz)] = (a->fieldZ[fidx(kx, ky, kz)]) * (c->field[fidx(kx, ky, kz)]);
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::add_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] + b->fieldX[fidx(kx, ky, kz)];
            fieldY[fidx(kx, ky, kz)] = a->fieldY[fidx(kx, ky, kz)] + b->fieldY[fidx(kx, ky, kz)];
            fieldZ[fidx(kx, ky, kz)] = a->fieldZ[fidx(kx, ky, kz)] + b->fieldZ[fidx(kx, ky, kz)];
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::sub_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = a->fieldX[fidx(kx, ky, kz)] - b->fieldX[fidx(kx, ky, kz)];
            fieldY[fidx(kx, ky, kz)] = a->fieldY[fidx(kx, ky, kz)] - b->fieldY[fidx(kx, ky, kz)];
            fieldZ[fidx(kx, ky, kz)] = a->fieldZ[fidx(kx, ky, kz)] - b->fieldZ[fidx(kx, ky, kz)];
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmVectorFieldOps::neg_plane(CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            fieldX[fidx(kx, ky, kz)] = -a->fieldX[fidx(kx, ky, kz)];
            fieldY[fidx(kx, ky, kz)] = -a->fieldY[fidx(kx, ky, kz)];
            fieldZ[fidx(kx, ky, kz)] = -a->fieldZ[fidx(kx, ky, kz)];
        }

    return 0;
}
