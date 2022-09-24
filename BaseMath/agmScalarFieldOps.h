#pragma once

class CagmVectorFieldOps;

//------------------------------------------------------------------
class CagmScalarFieldOps
{
friend class CagmVectorFieldOps;

public:
    int N[3];
    REALTYPE__ step[3];
    int * GetDimensions()
    {
        return N;
    }
    static uint32_t GetAllocSize(int *N)
    {
        return sizeof(REALTYPE__)*N[1]*N[2] + sizeof(CagmScalarFieldOps);
    }

protected:
    REALTYPE__ **field;
    int NphysL[3], NphysH[3];

public:
	CagmScalarFieldOps(int *_N, int *_DphysL = nullptr, int *_DphysH = nullptr);
	CagmScalarFieldOps(const CagmScalarFieldOps&);
	virtual ~CagmScalarFieldOps();

    CagmScalarFieldOps& operator=(const CagmScalarFieldOps&);

	uint32_t setDPhys(int *_DphysL, int *_DphysH);

    //REALTYPE__ getElement(int kx, int ky, int kz);
    REALTYPE__ *getAddress(int kx, int ky, int kz);

    virtual uint32_t SetMargins(CagmScalarFieldOps *source, int *Mmin, int *_DphysL = nullptr, int *_DphysH = nullptr);

    uint32_t SetSteps(REALTYPE__ *_step);
    REALTYPE__ *GetSteps();

	//uint32_t set(CagmScalarFieldOps *a);
	REALTYPE__ derivative(int kx, int ky, int kz, int dir);
    uint32_t div(CagmVectorFieldOps *a);
    uint32_t div31(CagmVectorFieldOps *a);
    uint32_t div42(CagmVectorFieldOps *a);
    uint32_t div41(CagmVectorFieldOps *a);
    uint32_t div5(CagmVectorFieldOps *a);
    uint32_t divScheme(CagmVectorFieldOps *a, int scheme);
    uint32_t dot(CagmVectorFieldOps *a, CagmVectorFieldOps *b, CagmVectorFieldOps *Weight = nullptr);
	uint32_t abs2(CagmVectorFieldOps *a, CagmVectorFieldOps *Weight = nullptr);
	uint32_t abs(CagmVectorFieldOps *a);
	uint32_t projection(CagmVectorFieldOps *a, REALTYPE__ *d);
	uint32_t inv(CagmScalarFieldOps *a);
	uint32_t inv(void);
	uint32_t mult(REALTYPE__ c, CagmScalarFieldOps *a);
	uint32_t mult(REALTYPE__ c);
	uint32_t mult(CagmScalarFieldOps *c, CagmScalarFieldOps *a);
	uint32_t mult(CagmScalarFieldOps *c);
    uint32_t add(CagmScalarFieldOps *a, CagmScalarFieldOps *b);
    uint32_t addPhys(CagmScalarFieldOps *a, CagmScalarFieldOps *b);
    uint32_t add(CagmScalarFieldOps *a);
    uint32_t addPhys(CagmScalarFieldOps *a);
    uint32_t sub(CagmScalarFieldOps *a, CagmScalarFieldOps *b);
    uint32_t sub(CagmScalarFieldOps *a);
    uint32_t neg(CagmScalarFieldOps *a);
    uint32_t neg();
    uint32_t acos();
    uint32_t power(CagmScalarFieldOps *a, REALTYPE__ pw);
    uint32_t power(REALTYPE__ pw);
    uint32_t zero();
    uint32_t zeroZ0();
    uint32_t setZlevel(int level, REALTYPE__ w);
    uint32_t setPlane(CagmScalarFieldOps *plane, int wplane, int from, int to);

	uint32_t LOS(CagmVectorFieldOps *a, REALTYPE__ *dircos);
    static uint32_t rotate2D(CagmScalarFieldOps *ac, CagmScalarFieldOps *as, CagmScalarFieldOps *acn, CagmScalarFieldOps *asn, REALTYPE__ cosz);
    uint32_t sqDiff(CagmScalarFieldOps *a1, CagmScalarFieldOps *a2);
    uint32_t relax(CagmScalarFieldOps *mult, CagmScalarFieldOps *weight);

    REALTYPE__ sum(CagmScalarFieldOps *weight = nullptr);
    REALTYPE__ sumPhys(CagmScalarFieldOps *weight = nullptr);
    REALTYPE__ sumPhysW(CagmScalarFieldOps *weight);
    REALTYPE__ avPhys(CagmScalarFieldOps *weight = nullptr);
    REALTYPE__ maxval(void);
    uint32_t limWeight(int limType, CagmScalarFieldOps *calc, CagmScalarFieldOps *cond);

protected:
	uint32_t Initialize(int *_N, int *_NphysL = nullptr, int *_NphysH = nullptr, REALTYPE__ *_step = nullptr);
  	uint32_t setNPhys(int *_NphysL, int *_NphysH);
    uint32_t Delete();
};

// static void indices(int p, int Np, int q, int Nq, int *p1, int *p2, REALTYPE__ *pf, int *q1, int *q2, REALTYPE__ *qf);

