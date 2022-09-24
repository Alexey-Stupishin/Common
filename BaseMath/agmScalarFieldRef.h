#if !defined(AGS_20150923_0653_AGMSCALARFIELDREF_H)
#define AGS_20150923_0653_AGMSCALARFIELDREF_H

#include "precompile.h"
#include "agmScalarFieldOps.h"

//------------------------------------------------------------------
class CagmScalarFieldRef : public CagmScalarFieldOps
{
public:

protected:

public:
	CagmScalarFieldRef(int *_N, CagmScalarFieldOps *source, int *M, int *Mmin)
        {
            init(_N);

            int ky, kz;
            for (ky = 0; ky < N[1]; ky++)
                for (kz = 0; kz < N[2]; kz++)
                    field[ky + kz*N[1]] = source->getAddress(Mmin[0], ky+Mmin[1], kz+Mmin[2]);
		}

	virtual ~CagmScalarFieldRef()
    {
    }

};

#endif AGS_20150923_0653_AGMSCALARFIELDREF_H
