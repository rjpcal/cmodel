///////////////////////////////////////////////////////////////////////
//
// doSimplex.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Sun Apr  1 19:52:50 2001
// written: Tue Feb 19 15:02:12 2002
// $Id$
//
//
// MATLAB Compiler: 2.1
// Date: Sun Apr  1 19:52:50 2001
//
// Arguments: "-B" "macro_default" "-O" "all" "-O"
// "fold_scalar_mxarrays:on" "-O" "fold_non_scalar_mxarrays:on" "-O"
// "optimize_integer_for_loops:on" "-O" "array_indexing:on" "-O"
// "optimize_conditionals:on" "-x" "-W" "mex" "-L" "C" "-t" "-T"
// "link:mexlibrary" "libmatlbmx.mlib" "-h" "doSimplex"
//
///////////////////////////////////////////////////////////////////////

#ifndef DOSIMPLEX_CC_DEFINED
#define DOSIMPLEX_CC_DEFINED

#include "mexbuf.h"
#include "mtx.h"
#include "mx.h"
#include "matlabfunction.h"
#include "simplexoptimizer.h"

#include "util/error.h"
#include "util/strings.h"

#include <iostream>
#include <libmatlbm.h>

#define LOCAL_PROF
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
          mexErrMsgTxt("Option must be an integer value "
                       "if not the default.");
        }
    }

  return Mx::getInt(arr);
}

namespace
{
  MexBuf* mexBuf = 0;
}

void InitializeModule_doSimplex(void)
{
  mexPrintf("loading 'doSimplex' mex file\n");

  mexBuf = new MexBuf;
#ifdef MIPS_PRO
  std::cout = mexBuf;
  std::cerr = mexBuf;
#else
  std::cout.rdbuf(mexBuf);
  std::cerr.rdbuf(mexBuf);
#endif
}

void TerminateModule_doSimplex(void)
{
  mexPrintf("unloading 'doSimplex' mex file\n");

  Util::Prof::printAtExit(false);

  // For some reason we get crashses if we try to reset these to their
  // original streambuf*'s; setting them to 0 effectively just shuts the
  // streams down
  std::cout.rdbuf(0);
  std::cerr.rdbuf(0);

  delete mexBuf;
}


/*
 * The function "MdoSimplex" is the implementation version of the "doSimplex"
 * M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/doSimplex.m" (lines 1-237).
 * It contains the actual compiled code for that M-function. It is a static
 * function and must only be called from one of the interface functions,
 * appearing below.
 */

static mxArray * MdoSimplex(mxArray * * fval,
                            mxArray * * exitflag_mx,
                            mxArray * * output,
                            mxArray * funfcn_mx,
                            mxArray * x_in,
                            mxArray * printtype_mx,
                            mxArray * tolx_mx,
                            mxArray * tolf_mx,
                            mxArray * maxfun_mx,
                            mxArray * maxiter_mx,
                            mxArray * debugFlags_mx,
                            int nvararg,
                            mxArray** pvararg)
{
DOTRACE("MdoSimplex");

  try
    {
#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
      if (Mx::hasField(debugFlags_mx, "debugFlag"))
        {
          int debugFlag = Mx::getIntField(debugFlags_mx, "debugFlag");

          if (debugFlag == -1)
            {
              Util::Prof::printAllProfData(std::cerr);
            }
          else if (debugFlag == -2)
            {
              Util::Prof::resetAllProfData();
            }

          return mxCreateScalarDouble(debugFlag);
        }
#endif

      const Mtx x(x_in, Mtx::BORROW);
      const int numModelParams = x.nelems();

      MatlabFunction objective(Mx::getString(funfcn_mx),
                               nvararg,
                               pvararg,
                               false);

      SimplexOptimizer opt(objective,
                           Mtx(x_in),
                           Mx::getString(printtype_mx),
                           numModelParams,
                           extractMaxIters(maxfun_mx, numModelParams),
                           extractMaxIters(maxiter_mx, numModelParams),
                           Mx::getDouble(tolx_mx),
                           Mx::getDouble(tolf_mx)
                           );

      int exitFlag = opt.optimize();

      *fval = mxCreateScalarDouble(opt.bestFval());
      *exitflag_mx = mxCreateScalarDouble(exitFlag);

      mlfIndexAssign(output, ".iterations",
                     mxCreateScalarDouble(opt.iterCount()));

      mlfIndexAssign(output, ".funcCount",
                     mxCreateScalarDouble(opt.funcCount()));

      mlfIndexAssign(output, ".algorithm", mxCreateString(opt.algorithm()));

      return opt.bestParams().makeMxArray();
    }
  catch (Util::Error& err)
    {
      mexErrMsgTxt(err.msg_cstr());
    }
  catch (...)
    {
      mexErrMsgTxt("an unknown C++ exception occurred.");
    }

  return (mxArray*) 0; // can't happen, but placate compiler
}

/*
 * The function "mlxDoSimplex" contains the feval interface for the "doSimplex"
 * M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/doSimplex.m" (lines 1-237).
 * The feval function calls the implementation version of doSimplex through
 * this function. This function processes any input arguments and passes them
 * to the implementation version of the function, appearing above.
 */
void mlxDoSimplex(int nlhs, mxArray * plhs[], int nrhs, mxArray * prhs[])
{
  const int NDECLARED = 8;

  mxArray* mprhs[NDECLARED];
  mxArray* mplhs[4];

  if (nlhs > 4)
    {
      mexErrMsgTxt("Error: doSimplex was called with more "
                   "than the declared number of outputs (4).");
    }
  if (nrhs < NDECLARED)
    {
      mexErrMsgTxt("Error: doSimplex was called with fewer "
                   "than the decalred number of inputs.");
    }

  for (int i = 0; i < 4; ++i)
    {
      mplhs[i] = 0;
    }
  for (int i = 0; i < NDECLARED; ++i)
    {
      mprhs[i] = prhs[i];
    }

  mlfEnterNewContext(0, 8,
                     mprhs[0], mprhs[1], mprhs[2], mprhs[3],
                     mprhs[4], mprhs[5], mprhs[6], mprhs[7]);

  int nvararg = nrhs - NDECLARED;
  mxArray** pvararg = prhs + NDECLARED;

  mplhs[0] = MdoSimplex(&mplhs[1],
                        &mplhs[2],
                        &mplhs[3],
                        mprhs[0],
                        mprhs[1],
                        mprhs[2],
                        mprhs[3],
                        mprhs[4],
                        mprhs[5],
                        mprhs[6],
                        mprhs[7],
                        nvararg,
                        pvararg);

  mlfRestorePreviousContext(0, 8,
                            mprhs[0], mprhs[1], mprhs[2], mprhs[3],
                            mprhs[4], mprhs[5], mprhs[6], mprhs[7]);

  plhs[0] = mplhs[0];
  for (int i = 1; i < 4 && i < nlhs; ++i)
    {
      if (i < nlhs)
        plhs[i] = mplhs[i];
      else
        mxDestroyArray(mplhs[i]);
    }
}

#ifndef MLF_V2
#define MLF_V2 1
#endif

static mexFunctionTableEntry function_table[1]
  = { { "doSimplex", mlxDoSimplex, -9, 4, (_mexLocalFunctionTable*)0 } };

static _mexInitTermTableEntry init_term_table[1]
  = { { InitializeModule_doSimplex, TerminateModule_doSimplex } };

static _mex_information _mex_info
  = { 1, 1, function_table, 0, NULL, 0, NULL, 1, init_term_table };

/*
 * The function "mexLibrary" is a Compiler-generated mex wrapper, suitable for
 * building a MEX-function. It initializes any persistent variables as well as
 * a function table for use by the feval function. It then calls the function
 * "mlxDoSimplex". Finally, it clears the feval table and exits.
 */
extern "C"
mex_information mexLibrary(void)
{
  return &_mex_info;
}

static const char vcid_doSimplex_cc[] = "$Header$";
#endif // !DOSIMPLEX_CC_DEFINED
