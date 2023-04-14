#include "stdDefinitions.h"
#include "agmLanczos.h"

//---------------------------------------------------------------------
CagmLanczos::CagmLanczos(int _a, int _hw)
        : a(_a)
        , hw(_hw)
    {
        L = a*hw + 1;
        halfKernel = new double[L];

        halfKernel[0] = 1;
        for (int i = 1; i < L; i++)
        {
            int x = (double)i/L * M_PI;
            halfKernel[i] = a * sin(x) * sin(x/a) /(x*x);
        }

//                    this.lanz_step = 1/this.N;

    }

//---------------------------------------------------------------------
unsigned int CagmLanczos::resample(int N, double *from, int Nto, double *to)
{
    double step = (double)(Nto-1)/(double)(N-1);
    double *xx = new double[Nto];

    memset(to, 0, Nto*sizeof(double));
    for (int i = 0; i < L; i++)
    {
        double x = i*step;
        for (int aa = -a; aa <= a; aa++)
        {
            int idxf = floor(x) + aa;
            if (idxf < 0 || idxf >= N)
                continue;
            int idxl = floor(fabs(x - idxf)*L);
            if (idxl < 0 || idxf >= L)
                continue;
            to[i] += from[idxf]*halfKernel[idxl];
        }
    }

    delete[] xx;

    return 0;
}

//---------------------------------------------------------------------
unsigned int CagmLanczos::resample(int dim, int *N, double *from, int *Nto, double *to)
{
    int *c_Nfrom = new int[dim];
    double *c_from ;

    int *c_Nto = new int[dim];
    double *c_to;

    // preface
    c_from = from;
    memcpy(c_Nfrom, N, dim*sizeof(int));
    memcpy(c_Nto, Nto, dim*sizeof(int));

    for (int d = 0; d < dim; d++)
    {
        // set
        c_Nto[d] = Nto[d];
        if (d == dim-1)
            c_to = to;
        else
            c_to = new double[];

        // process

        // finish
        if (d < dim - 1)
            delete[] c_to;

    }

    return 0;
}
