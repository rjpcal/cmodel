///////////////////////////////////////////////////////////////////////
//
// annealHelper.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Tue Feb 19 10:20:05 2002
// $Id$
//
//
// MATLAB Compiler: 2.1
// Date: Fri Mar 23 17:17:00 2001
//
// Arguments: "-B" "macro_default" "-O" "all" "-O" "fold_scalar_mxarrays:on"
// "-O" "fold_non_scalar_mxarrays:on" "-O" "optimize_integer_for_loops:on"
// "-O" "array_indexing:on" "-O" "optimize_conditionals:on" "-x" "-W" "mex"
// "-L" "C" "-t" "-T" "link:mexlibrary" "libmatlbmx.mlib" "-h"
// "annealHelper"
//
///////////////////////////////////////////////////////////////////////

#ifndef ANNEALHELPER_CC_DEFINED
#define ANNEALHELPER_CC_DEFINED

#include "annealHelper.h"

#include "annealingoptimizer.h"
#include "matlabfunction.h"
#include "mexbuf.h"
#include "mx.h"

#include "util/error.h"

#include <iostream>
#include <libmatlbm.h>

#include "util/trace.h"

namespace
{
  MexBuf* mexBuf = 0;

  std::streambuf* coutOrigBuf = 0;
  std::streambuf* cerrOrigBuf = 0;

  void debugDestroyArray(mxArray* arr, const char* name,
                         const char* file, int line)
  {
    mexPrintf("file %s, line %d, about to destroy %s...", file, line, name);
    mxDestroyArray(arr);
    mexPrintf(" done.\n");
  }
}

#ifdef MXDEBUG
#define mxDestroyArray(x)  debugDestroyArray(x, #x, __FILE__, __LINE__);
#endif

void InitializeModule_annealHelper()
{
  mexPrintf("loading 'annealHelper' mex file\n");

  mexBuf = new MexBuf;
#ifdef MIPS_PRO
  std::cout = mexBuf;
  std::cerr = mexBuf;
#else
  coutOrigBuf = std::cout.rdbuf(mexBuf);
  cerrOrigBuf = std::cerr.rdbuf(mexBuf);
#endif
}

void TerminateModule_annealHelper()
{
  mexPrintf("unloading 'annealHelper' mex file\n");

  Util::Prof::printAtExit(false);

  // For some reason we get crashses if we try to reset these to their
  // original streambuf*'s; setting them to 0 effectively just shuts the
  // streams down
  std::cout.rdbuf(0);
  std::cerr.rdbuf(0);

  delete mexBuf;
}

/*
 * The function "mlxAnnealHelper" contains the feval interface for the
 * "annealHelper" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealHelper.m" (lines
 * 1-28). The feval function calls the implementation version of annealHelper
 * through this function. This function processes any input arguments and
 * passes them to the implementation version of the function, appearing above.
 */
void mlxAnnealHelper(int nlhs, mxArray* plhs[], int nrhs, mxArray* prhs[])
{
DOTRACE("mlxAnnealHelper");

  if (nlhs > 1)
    {
      mexErrMsgTxt("Error: annealHelper was called with more "
                   "than the declared number of outputs (1).");
    }

  const int NDECLARED = 2;

  if (nrhs < NDECLARED)
    {
      mexErrMsgTxt("Error: annealHelper was called with fewer "
                   "than the declared number of inputs (7).");
    }

  mlfEnterNewContext(0, NDECLARED, prhs[0], prhs[1]);

  int nvararg = nrhs - NDECLARED;
  mxArray** pvararg = prhs + NDECLARED;

  try
    {
      AnnealOpts opts(prhs[0]);

      MatlabFunction objective(Mx::getString(prhs[1]), // funcName
                               nvararg,
                               pvararg,
                               opts.canUseMatrix);

      AnnealingOptimizer ar(objective, opts);

      ar.optimize();

      plhs[0] = ar.getOutput();
    }
  catch (Util::Error& err)
    {
      mexErrMsgTxt(err.msg_cstr());
    }
  catch (...)
    {
      mexErrMsgTxt("an unknown C++ exception occurred.");
    }

  mlfRestorePreviousContext(0, NDECLARED, prhs[0], prhs[1]);
}

static const char vcid_annealHelper_cc[] = "$Header$";
#endif // !ANNEALHELPER_CC_DEFINED
