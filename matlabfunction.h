///////////////////////////////////////////////////////////////////////
//
// matlabfunction.h
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 18 11:07:28 2002
// written: Tue Feb 19 09:47:00 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MATLABFUNCTION_H_DEFINED
#define MATLABFUNCTION_H_DEFINED

#include "multivarfunction.h"

#include "mtx.h"

#include "util/error.h"
#include "util/strings.h"

#include <libmatlb.h>

///////////////////////////////////////////////////////////////////////
/**
 *
 * MatlabFunction represents a function that takes one variable argument, plus
 * any number of trailing arguments that are bound at the time the object is
 * created.
 *
 **/
///////////////////////////////////////////////////////////////////////

class Objective : public MultivarFunction
{
private:
  const fstring itsFuncName;
  const int itsNvararg;
  mxArray** const itsPvararg;
  const bool itsCanUseMatrix;

  mxArray** itsPrhs;

public:
  Objective(const fstring& funcName, int nvararg, mxArray** pvararg,
            bool canUseMatrix = false)
    :
    itsFuncName(funcName),
    itsNvararg(nvararg),
    itsPvararg(pvararg),
    itsCanUseMatrix(canUseMatrix),
    itsPrhs(new (mxArray*)[nvararg+1])
  {
    itsPrhs[0] = 0;

    for (int i = 0; i < itsNvararg; ++i)
      {
        itsPrhs[i+1] = 0;
        mlfAssign(&itsPrhs[i+1], mxDuplicateArray(itsPvararg[i]));
      }
  }

  virtual ~Objective()
  {
    for (int i = 0; i < itsNvararg; ++i)
      {
        mxDestroyArray(itsPrhs[i+1]);
      }

    delete [] itsPrhs;
  }

private:
  Mtx evaluateParallel(const Mtx& models) const
  {
//     DOTRACE("evaluateParallel");

    mxArray* costs_mx = 0;

    // Since there is no mlfAssign here, this is treated as an unbound
    // (temporary) array, which is automatically destroyed during
    // mexCallMATLAB() below
    itsPrhs[0] = models.makeMxArray();

    int err =
      mexCallMATLAB(1, &costs_mx, itsNvararg+1, itsPrhs, itsFuncName.c_str());

    itsPrhs[0] = 0;

    if (err != 0)
      throw Util::Error("mexCallMATLAB failed in evaluateParallel");

    Mtx costs(costs_mx, Mtx::COPY);

    mxDestroyArray(costs_mx);

    return costs;
  }

  Mtx evaluateSerial(const Mtx& models) const
  {
//     DOTRACE("evaluateSerial");

    const int numModels = models.ncols();

    Mtx result(numModels, 1);

    mxArray* cost_mx = 0;

    for (int e = 0; e < numModels; ++e)
      {
        Mtx currentModel = models.column(e);

        // Since there is no mlfAssign here, this is treated as an unbound
        // (temporary) array, which is automatically destroyed during
        // mexCallMATLAB() below
        itsPrhs[0] = currentModel.makeMxArray();

        int err =
          mexCallMATLAB(1, &cost_mx, itsNvararg+1, itsPrhs, itsFuncName.c_str());

        itsPrhs[0] = 0;

        if (err != 0)
          throw Util::Error("mexCallMATLAB failed in evaluateSerial");

        result.at(e) = mxGetScalar(cost_mx);

        mxDestroyArray(cost_mx);
      }

    return result;
  }

protected:
  virtual Mtx doEvaluateEach(const Mtx& models)
  {
    if (itsCanUseMatrix)
      return evaluateParallel(models);

    return evaluateSerial(models);
  }

  virtual double doEvaluate(const Mtx& model)
  {
    Mtx result = doEvaluateEach(model.asColumn());
    return result.at(0);
  }
};

static const char vcid_matlabfunction_h[] = "$Header$";
#endif // !MATLABFUNCTION_H_DEFINED
