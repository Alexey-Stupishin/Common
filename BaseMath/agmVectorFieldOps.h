#pragma once

#include "math.h"
#include "agmRKF45.h"

class CagmScalarFieldOps;
class CagmRotate3D;

//------------------------------------------------------------------
class CagmVectorFieldOps
{
friend class CagmScalarFieldOps;
//friend CagmVectorFieldOps& operator*(const CagmVectorFieldOps&, const CagmVectorFieldOps&);

public:
    typedef enum {None = 0, Boundary = 2, OutOfCube = 4, BufferOverload = 8, RKF45Problem = 16
                 }  Status;
    typedef enum {Linear = 0, Lanczos = 1
                 }  Interpolator;

    static REALTYPE__ __norm(REALTYPE__ *v1, REALTYPE__ *v2)
    {
        return sqrt(v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
    }

public:
    int N[3];
    REALTYPE__ step[3];
    int * GetDimensions()
    {
        return N;
    }
    static uint32_t GetAllocSize(int *N)
    {
        return 3*sizeof(REALTYPE__)*N[1]*N[2] + sizeof(CagmVectorFieldOps);
    }

public:
    REALTYPE__ **fieldX, **fieldY, **fieldZ;
protected:
    int NphysL[3], NphysH[3];

public:
	CagmVectorFieldOps(int *_N, int *_DphysL = nullptr, int *_DphysH = nullptr);
	CagmVectorFieldOps(const CagmVectorFieldOps&);
	virtual ~CagmVectorFieldOps();

    CagmVectorFieldOps& operator=(const CagmVectorFieldOps&);

	uint32_t setDPhys(int *_DphysL, int *_DphysH);

    REALTYPE__ *getAddress(int v, int kx, int ky, int kz);

	virtual uint32_t SetMargins(CagmVectorFieldOps *source, int *Mmin, int *_DphysL = nullptr, int *_DphysH = nullptr);
    uint32_t setRefField(REALTYPE__ *X, REALTYPE__ *Y, REALTYPE__ *Z, int *M);

    uint32_t SetSteps(REALTYPE__ *_step);
    REALTYPE__ *GetSteps();

	uint32_t cross(CagmVectorFieldOps *a, const CagmVectorFieldOps *b);
    //CagmVectorFieldOps& operator*(const CagmVectorFieldOps&);
    uint32_t blockF(CagmVectorFieldOps *B, CagmScalarFieldOps *N, CagmScalarFieldOps *w, CagmVectorFieldOps *gradw);
    uint32_t cross(CagmVectorFieldOps *a);
	uint32_t rot(CagmVectorFieldOps *a);
    uint32_t rot31(CagmVectorFieldOps *a);
    uint32_t rot42(CagmVectorFieldOps *a);
    uint32_t rot41(CagmVectorFieldOps *a);
    uint32_t rotScheme(CagmVectorFieldOps *a, int scheme);
    uint32_t rot5(CagmVectorFieldOps *a);
	uint32_t grad(CagmScalarFieldOps *a);
    uint32_t grad31(CagmScalarFieldOps *a);
    uint32_t grad42(CagmScalarFieldOps *a);
    uint32_t grad41(CagmScalarFieldOps *a);
    uint32_t gradScheme(CagmScalarFieldOps *a, int scheme);
    uint32_t grad5(CagmScalarFieldOps *a);
	uint32_t mult(REALTYPE__ c, CagmVectorFieldOps *a);
	uint32_t mult(REALTYPE__ c);
	uint32_t mult(CagmScalarFieldOps *c, CagmVectorFieldOps *a);
	uint32_t multPhys(CagmScalarFieldOps *c, CagmVectorFieldOps *a);
	uint32_t multPhys(CagmVectorFieldOps *c, CagmVectorFieldOps *a);
	uint32_t mult(CagmVectorFieldOps *a, CagmScalarFieldOps *c);
	uint32_t multPhys(CagmVectorFieldOps *a, CagmScalarFieldOps *c);
	uint32_t mult(CagmScalarFieldOps *c);
	uint32_t multPhys(CagmScalarFieldOps *c);
	uint32_t multPhys(CagmVectorFieldOps *c);
	uint32_t mult(CagmVectorFieldOps *a, REALTYPE__ *d);
    uint32_t add(CagmVectorFieldOps *a, CagmVectorFieldOps *b);
    uint32_t addPhys(CagmVectorFieldOps *a, CagmVectorFieldOps *b);
    uint32_t add(CagmVectorFieldOps *a);
    uint32_t addPhys(CagmVectorFieldOps *a);
	uint32_t sub(CagmVectorFieldOps *a, CagmVectorFieldOps *b);
	uint32_t subPhys(CagmVectorFieldOps *a, CagmVectorFieldOps *b);
    uint32_t sub(CagmVectorFieldOps *a);
    uint32_t subPhys(CagmVectorFieldOps *a);
    uint32_t neg(CagmVectorFieldOps *a);
    uint32_t neg();
    uint32_t zero();

    uint32_t getThetaMetrics(REALTYPE__ relBound, REALTYPE__ *theta, REALTYPE__ *thetaj);

    //WiegelmannWeightBound

    uint32_t shift(int n);
    uint32_t setZlevel(int wplane, int level, REALTYPE__ w = 1.0);
    uint32_t setVector(REALTYPE__ *w);

    uint32_t getPlane(CagmVectorFieldOps *plane, int wplane, int from, int to);
    uint32_t setPlaneComp(CagmVectorFieldOps *plane, int wplane, int wcomp, int from, int to);
    uint32_t setPlane(CagmVectorFieldOps *plane, int wplane, int from, int to);

    uint32_t getComponent(CagmScalarFieldOps *comp, int wcomp);
    uint32_t setComponent(CagmScalarFieldOps *comp, int wcomp);
    uint32_t getTransv(CagmScalarFieldOps *a);

    uint32_t getBounds(CagmVectorFieldOps *boundsx, CagmVectorFieldOps *boundsy, CagmVectorFieldOps *boundsz);
    uint32_t setBounds(CagmVectorFieldOps *boundsx, CagmVectorFieldOps *boundsy, CagmVectorFieldOps *boundsz);

	uint32_t stretch(CagmVectorFieldOps *a, Interpolator = Linear, REALTYPE__ p1 = 0, REALTYPE__ p2 = 0, REALTYPE__ p3 = 0);
    uint32_t conv(CagmVectorFieldOps *src, CagmScalarFieldOps *win);
    uint32_t inCube(const REALTYPE__ * coord, const REALTYPE__ absBoundAchieve = 0, const REALTYPE__ relBoundAchieve = 0);
    uint32_t getPoint(const REALTYPE__ *coord, REALTYPE__ *vect);
    //Status getLine(REALTYPE__ *start, int ninpt, int dir, REALTYPE__ step, REALTYPE__ tolerance, REALTYPE__ boundAchieve, 
    //               int maxResult, int *lengthLines, REALTYPE__ *coord, int *codes);
    //Status getLine1(REALTYPE__ *start, int direction, int maxLength, REALTYPE__ step, REALTYPE__ tolerance, REALTYPE__ boundAchieve, 
    //                int& length, REALTYPE__ *coord, REALTYPE__ *steps);

    Status getOneFullLine(CagmRKF45 *rkf45, REALTYPE__ *start, int direction, REALTYPE__ step, REALTYPE__ boundAchieve, REALTYPE__ boundAchieveBottom,
        int maxResult, int *length, REALTYPE__ *coord, int *status);

    uint32_t rotate3D(CagmRotate3D *, bool);
    uint32_t planeDerivative(int layer, REALTYPE__ *d);
    uint32_t planeDerivative2(int layer, REALTYPE__ *d, REALTYPE__ *dd);

protected:
	uint32_t Initialize(int *_N, int *_NphysL = nullptr, int *_NphysH = nullptr, REALTYPE__ *_step = nullptr);
  	uint32_t setNPhys(int *_NphysL, int *_NphysH);
    Status getOneLine(CagmRKF45 *rkf45, CagmRKF45Vect *rkfv, REALTYPE__ step, REALTYPE__ *coord, int maxlen, int *length, CagmRKF45::Status *status, bool noDuplicate = false);
    uint32_t Delete();
};
