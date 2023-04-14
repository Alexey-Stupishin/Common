#pragma once

#define _USE_MATH_DEFINES
#include <math.h>

class CagmLanczos
{
protected:
    int a;
    int hw;
    int L;
    double *halfKernel;

public:
    CagmLanczos(int _a = 2, int _hw = 1000000);
    virtual ~CagmLanczos()
    {
        delete [] halfKernel;
    }

    unsigned int resample(int N, double *from, int Nto, double *to);
    unsigned int resample(int dim, int *N, double *from, int *Nto, double *to);
};

