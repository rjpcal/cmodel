///////////////////////////////////////////////////////////////////////
//
// annealHelper.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Tue Sep 28 13:17:40 2004
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef ANNEALHELPER_CC_DEFINED
#define ANNEALHELPER_CC_DEFINED

#include "cmodel/annealingoptimizer.h"
#include "cmodel/matlabfunction.h"

#include "mx/mexpkg.h"
#include "mx/mx.h"

#include "util/error.h"

#include "util/trace.h"

void annealHelper(int nlhs, mxArray* plhs[],
                  int nrhs, const mxArray* prhs[])
{
DOTRACE("annealHelper");

  const int NDECLARED = 2; // The minimum number of input arguments

  int nvararg = nrhs - NDECLARED;
  const mxArray** pvararg = prhs + NDECLARED;

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
  if (nvararg > 0 && Mx::hasField(pvararg[0], "debugFlag"))
    {
      int debugFlag = Mx::getIntField(pvararg[0], "debugFlag");

      if (debugFlag == -1)       Util::Prof::printAllProfData(std::cerr);
      else if (debugFlag == -2)  Util::Prof::resetAllProfData();

      plhs[0] = mxCreateScalarDouble(debugFlag);
    }
  else
    {
#endif
      AnnealOpts opts(prhs[0]);

      MatlabFunction objective(Mx::getString(prhs[1]), // funcName
                               nvararg,
                               pvararg,
                               opts.canUseMatrix);

      AnnealingOptimizer ar(objective, opts);

      ar.optimize();

      plhs[0] = ar.getOutput();
#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
    }
#endif
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

static const char vcid_annealHelper_cc[] = "$Header$";
#endif // !ANNEALHELPER_CC_DEFINED
