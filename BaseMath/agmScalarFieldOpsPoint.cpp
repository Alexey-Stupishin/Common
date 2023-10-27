#include "stdDefinitions.h"
#include <math.h>

#include "agsFieldsCommon.h"
#include "agmScalarFieldOps.h"
#include "agmVectorFieldOps.h"

#define fidx(kx, ky, kz) (ky)+(kz)*N[1]][(kx)

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::div_plane(CagmVectorFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    REALTYPE_A dx, dy, dz;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            if (kx == 0)
                dx = -3 * a->fieldX[fidx(0, ky, kz)] + 4 * a->fieldX[fidx(1, ky, kz)] - a->fieldX[fidx(2, ky, kz)];
            else if (kx == N[0] - 1)
                dx = a->fieldX[fidx(N[0] - 3, ky, kz)] - 4 * a->fieldX[fidx(N[0] - 2, ky, kz)] + 3 * a->fieldX[fidx(N[0] - 1, ky, kz)];
            else
                dx = a->fieldX[fidx(kx + 1, ky, kz)] - a->fieldX[fidx(kx - 1, ky, kz)];

            if (ky == 0)
                dy = -3 * a->fieldY[fidx(kx, 0, kz)] + 4 * a->fieldY[fidx(kx, 1, kz)] - a->fieldY[fidx(kx, 2, kz)];
            else if (ky == N[1] - 1)
                dy = a->fieldY[fidx(kx, N[1] - 3, kz)] - 4 * a->fieldY[fidx(kx, N[1] - 2, kz)] + 3 * a->fieldY[fidx(kx, N[1] - 1, kz)];
            else
                dy = a->fieldY[fidx(kx, ky + 1, kz)] - a->fieldY[fidx(kx, ky - 1, kz)];

            if (kz == 0)
                dz = -3 * a->fieldZ[fidx(kx, ky, 0)] + 4 * a->fieldZ[fidx(kx, ky, 1)] - a->fieldZ[fidx(kx, ky, 2)];
            else if (kz == N[2] - 1)
                dz = a->fieldZ[fidx(kx, ky, N[2] - 3)] - 4 * a->fieldZ[fidx(kx, ky, N[2] - 2)] + 3 * a->fieldZ[fidx(kx, ky, N[2] - 1)];
            else
                dz = a->fieldZ[fidx(kx, ky, kz + 1)] - a->fieldZ[fidx(kx, ky, kz - 1)];

            field[fidx(kx, ky, kz)] = (dx*step[0] + dy*step[1] + dz*step[2]) * 0.5;
            //field[fidx(kx, ky, kz)] = (dx + dy + dz) * 0.5;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::dot_plane(CagmVectorFieldOps *a, CagmVectorFieldOps *b, int kz, CagmScalarFieldOps *Weight)
{
    // check equiv. sizes!
    REALTYPE_A w = 1;
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            if (Weight)
                w = Weight->field[fidx(kx, ky, kz)];
            field[fidx(kx, ky, kz)] = (a->fieldX[fidx(kx, ky, kz)] * b->fieldX[fidx(kx, ky, kz)]
                                     + a->fieldY[fidx(kx, ky, kz)] * b->fieldY[fidx(kx, ky, kz)]
                                     + a->fieldZ[fidx(kx, ky, kz)] * b->fieldZ[fidx(kx, ky, kz)]
                ) * w;
        }

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::abs2_plane(CagmVectorFieldOps *a, int kz, CagmScalarFieldOps *Weight)
{
    return dot_plane(a, a, kz, Weight);
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::inv_plane(CagmScalarFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
            if (a->field[fidx(kx, ky, kz)] != 0)
                field[fidx(kx, ky, kz)] = 1.0 / a->field[fidx(kx, ky, kz)];
            else
                field[fidx(kx, ky, kz)] = 0;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult_plane(REALTYPE_A c, CagmScalarFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
            field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] * c;

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::mult_plane(CagmScalarFieldOps *c, CagmScalarFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
            if (a->field[fidx(kx, ky, kz)] == 0)
                field[fidx(kx, ky, kz)] = 0;
            else
                field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] * c->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::add_plane(CagmScalarFieldOps *a, CagmScalarFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
            field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] + b->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::sub_plane(CagmScalarFieldOps *a, CagmScalarFieldOps *b, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
            field[fidx(kx, ky, kz)] = a->field[fidx(kx, ky, kz)] - b->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
uint32_t CagmScalarFieldOps::neg_plane(CagmScalarFieldOps *a, int kz)
{
    // check equiv. sizes!
    int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
            field[fidx(kx, ky, kz)] = -a->field[fidx(kx, ky, kz)];

    return 0;
}

//-----------------------------------------------------------------------
REALTYPE_A CagmScalarFieldOps::sum_plane(int kz, CagmScalarFieldOps *weight)
{
    REALTYPE_A wsum = 0;
	int kx, ky;
    for (ky = 0; ky < N[1]; ky++)
        for (kx = 0; kx < N[0]; kx++)
        {
            if (weight)
    			wsum += field[fidx(kx, ky, kz)] * weight->field[fidx(kx, ky, kz)];
            else
    			wsum += field[fidx(kx, ky, kz)];
        }

    return wsum;
}
