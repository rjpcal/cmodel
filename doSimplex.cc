///////////////////////////////////////////////////////////////////////
//
// doSimplex.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Sun Apr  1 19:52:50 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef DOSIMPLEX_CC_DEFINED
#define DOSIMPLEX_CC_DEFINED

#include "mtx/mtx.h"

#include "mx/mexpkg.h"
#include "mx/mx.h"
#include "mx/mxwrapper.h"

#include "optim/matlabfunction.h"
#include "optim/simplexoptimizer.h"

#include "util/error.h"

#include "util/trace.h"

int extractMaxIters(const mxArray* arr, int numModelParams)
{
  if (mxIsChar(arr))
    {
      if (Mx::getString(arr) == "200*numberofvariables")
        {
          return 200*numModelParams;
        }
      else
        {
          throw Util::Error("Option must be an integer value "
                            "if not the default.");
        }
    }

  return Mx::getInt(arr);
}

void doSimplex(int nlhs, mxArray* plhs[],
               int nrhs, const mxArray* prhs[])
{
  const mxArray* funfcn_mx     =  prhs[0];
  const mxArray* x_in          =  prhs[1];
  const mxArray* printtype_mx  =  prhs[2];
  const mxArray* tolx_mx       =  prhs[3];
  const mxArray* tolf_mx       =  prhs[4];
  const mxArray* maxfun_mx     =  prhs[5];
  const mxArray* maxiter_mx    =  prhs[6];
  const mxArray* debugFlags_mx =  prhs[7];

  (void) debugFlags_mx; // This variable is no longer used, but we
                        // need to make sure that any matlab scripts
                        // that use this mex function are updated
                        // before we remove 'debugFlags' from the
                        // function prototype

  const int NDECLARED = 8;

  int nvararg = nrhs - NDECLARED;
  const mxArray** pvararg = prhs + NDECLARED;

  const mtx x(x_in, mtx::BORROW);
  const int numModelParams = x.nelems();

  // Setup the objective function
  MatlabFunction objective(Mx::getString(funfcn_mx),
                           nvararg,
                           pvararg,
                           false);

  // Create the optimizer...
  SimplexOptimizer opt(objective,
                       mtx(x_in),
                       Mx::getString(printtype_mx),
                       numModelParams,
                       extractMaxIters(maxfun_mx, numModelParams),
                       extractMaxIters(maxiter_mx, numModelParams),
                       Mx::getDouble(tolx_mx),
                       Mx::getDouble(tolf_mx)
                       );

  // ... and let it do its optimization
  const int exitFlag = opt.optimize();

  // Now return results back to the matlab caller

  plhs[0] = opt.bestParams().makeMxArray();

  plhs[1] = mxCreateScalarDouble(opt.bestFval());
  plhs[2] = mxCreateScalarDouble(exitFlag);

  MxWrapper output_wrapper(plhs[3]);

  output_wrapper.setField("iterations",
                          mxCreateScalarDouble(opt.iterCount()));

  output_wrapper.setField("funcCount",
                          mxCreateScalarDouble(opt.funcCount()));

  output_wrapper.setField("algorithm",
                          mxCreateString(opt.algorithm()));
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
      mexPkg = new MexPkg("doSimplex", &terminateModule);
      mexPkg->addFcn("doSimplex", &doSimplex, -9, 4);
    }

  mexPkg->invokeFcn(nlhs, plhs, nrhs, prhs);
}

static const char vcid_doSimplex_cc[] = "$Header$";
#endif // !DOSIMPLEX_CC_DEFINED
