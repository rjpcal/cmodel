///////////////////////////////////////////////////////////////////////
//
// annealHelper.cc
//
// Copyright (c) 2001-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar 23 17:17:00 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef ANNEALHELPER_CC_DEFINED
#define ANNEALHELPER_CC_DEFINED

#include "mx/mexpkg.h"
#include "mx/mx.h"
#include "mx/mxwrapper.h"

#include "optim/annealingoptimizer.h"
#include "optim/matlabfunction.h"

#include "util/trace.h"

void annealHelper(int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[])
{
DOTRACE("annealHelper");

  const int NDECLARED = 2; // The minimum number of input arguments

  int nvararg = nrhs - NDECLARED;
  const mxArray** pvararg = prhs + NDECLARED;

  AnnealOpts opts((mx_wrapper(prhs[0])));

  MatlabFunction objective(mx::as_string(prhs[1]), // funcName
                           nvararg,
                           pvararg,
                           opts.canUseMatrix);

  AnnealingOptimizer ar(objective, opts);

  ar.optimize();

  mx_wrapper result;
  ar.getOutput(result);
  plhs[0] = result.release();
}

namespace
{
  MexPkg* mexPkg = 0;

  void terminateModule()
  {
  DOTRACE("terminateModule");
    delete mexPkg;
    mexPkg = 0;
  }
}

extern "C"
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
DOTRACE("mexFunction");

  if (mexPkg == 0)
    {
      mexPkg = new MexPkg("annealHelper", terminateModule);
      mexPkg->addFcn("annealHelper", &annealHelper, -2, 1);
    }

  mexPkg->invokeFcn(nlhs, plhs, nrhs, prhs);
}

static const char vcid_annealHelper_cc[] = "$Id$ $URL$";
#endif // !ANNEALHELPER_CC_DEFINED
