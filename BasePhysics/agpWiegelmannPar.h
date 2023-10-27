#pragma once

#include "stdDefinitions.h"
#include "TaskQueueProcessor.h"

class CagmScalarFieldOps;
class CagmVectorFieldOps;
class CagmScalarField;
class CagmVectorField;

//---------------------------------------------------------------------------------------
class PCOTask : public ATQPTask
{
public:
    int z;

public:
    PCOTask(int _id) : ATQPTask(_id) {}
    void setData(int _z)
    {
        z = _z;
    }
};

//---------------------------------------------------------------------------------------
class PCOTaskFactory : public ATQPTaskFactory
{
public:
    PCOTaskFactory() : ATQPTaskFactory() {}
    virtual ~PCOTaskFactory() {}

    virtual ATQPTask *create()
    {
        return new PCOTask(counter++);
    }
};

//---------------------------------------------------------------------------------------
class PCOSupervisor : public ATQPSupervisor
{
public:
    PCOSupervisor(int *N, ATQPTaskFactory* _factory, ATQPSynchonizer *_sync)
    : ATQPSupervisor(N[2], _factory, _sync)
    {
        int cnt = 0;
        for (int z = 0; z < N[2]; z++)
            ((PCOTask *)tasks[cnt++])->setData(z);
    }

    virtual ~PCOSupervisor() {}
};

class CagpWiegelmann;
//---------------------------------------------------------------------------------------
class PCOProcessor : public ATQPProcessor
{
protected:
    PCOTask *task;
    CagpWiegelmann* host;
    int stencil;

public:
    PCOProcessor(int id, ATQPSynchonizer *_sync, CagpWiegelmann* _host, int _stencil) 
        : ATQPProcessor(id, _sync)
        , host(_host)
        , stencil(_stencil)
    {}
    virtual ~PCOProcessor() {}

    virtual bool setTask(ATQPTask * _task);
    virtual bool proceed();
};

//------------------------------------------------------------------
class CagpWiegelmann // : public CParCubeOps
{
friend class PCOProcessor;

protected:
    int N[3];
    REALTYPE_A dircos[3];

    // sources input
    CagmVectorField *vB;
    CagmScalarField *sW;
    double *vdircos;

    // sources output
    CagmVectorFieldOps *vF;

    // sources constrains
    CagmVectorField *baseField;
    CagmVectorField *baseWeight;
    CagmScalarField *absField;
    CagmScalarField *absWeight;
    CagmScalarField *losField;
    CagmScalarField *losWeight;

    // calculate
    CagmVectorField *vgradW;

    // intermediate
  	CagmVectorField *rotB, *Wa, *Wb, *WaxB, *v;
	CagmScalarField *B2, *divB, *Wa2Wb2, *WbxB, *s;

    PCOSupervisor *supervisor;
    std::vector<ATQPProcessor *> processors;
    TaskQueueProcessor *main_proc;
    PCOTaskFactory *factory;
    int nProc;

    int mode;

    double *L;

    int iterN;
    int depth;

protected:

public:
    CagpWiegelmann(int *_N, int _n_threads, int _stencil
        , CagmVectorField *_sourceB, CagmScalarField *_sourceW
        , double *_vcos
        , CagmVectorFieldOps *_outF
        , CagmVectorField *_baseField, CagmVectorField *_baseWeight
        , CagmScalarField *_absField, CagmScalarField *_absWeight
        , CagmScalarField *_losField, CagmScalarField *_losWeight);
	virtual ~CagpWiegelmann();

    double step(int _iterN, int _depth);
};
