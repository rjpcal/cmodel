///////////////////////////////////////////////////////////////////////
//
// matlabfunction.cc
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 18 11:10:06 2002
// written: Tue Sep 28 13:14:28 2004
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MATLABFUNCTION_CC_DEFINED
#define MATLABFUNCTION_CC_DEFINED

#include "cmodel/matlabfunction.h"

#include "mtx/mtx.h"

#include "mx/mx.h"

#include "util/error.h"

#include <matrix.h>
#include <mex.h>

#include "util/trace.h"

MatlabFunction::MatlabFunction(const fstring& funcName,
                               int nvararg,
                               const mxArray** pvararg,
                               bool canUseMatrix)
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
      itsPrhs[i+1] = mxDuplicateArray(itsPvararg[i]);
    }
}

MatlabFunction::~MatlabFunction()
{
  for (int i = 0; i < itsNvararg; ++i)
    {
      mxDestroyArray(itsPrhs[i+1]);
    }

  delete [] itsPrhs;
}

Mtx MatlabFunction::evaluateParallel(const Mtx& models) const
{
DOTRACE("MatlabFunction::evaluateParallel");

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

Mtx MatlabFunction::evaluateSerial(const Mtx& models) const
{
DOTRACE("MatlabFunction::evaluateSerial");

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

      result.at(e) = Mx::getDouble(cost_mx);

      mxDestroyArray(cost_mx);
    }

  return result;
}

Mtx MatlabFunction::doEvaluateEach(const Mtx& models)
{
  if (itsCanUseMatrix)
    return evaluateParallel(models);

  return evaluateSerial(models);
}

double MatlabFunction::doEvaluate(const Mtx& model)
{
  Mtx result = doEvaluateEach(model.asColumn());
  return result.at(0);
}

static const char vcid_matlabfunction_cc[] = "$Header$";
#endif // !MATLABFUNCTION_CC_DEFINED
