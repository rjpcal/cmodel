///////////////////////////////////////////////////////////////////////
//
// doSimplex.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Sun Apr  1 19:52:50 2001
// written: Thu Feb 14 14:13:39 2002
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

#include "doSimplex.h"

#include "util/error.h"
#include "mexbuf.h"
#include "mtx.h"
#include "mxwrapper.h"
#include "multivarfunction.h"
#include "simplexoptimizer.h"
#include "util/strings.h"

#include <iostream.h>
#include <iomanip.h>
#include "libmatlbm.h"

#define LOCAL_PROF
#include "util/trace.h"

class MatlabFunction : public MultivarFunction
{
  mxArray* itsFunfcn;
  mxArray* itsVarargin_ref;

  static mxArray* getref(mxArray* varargin)
    {
      return mlfIndexRef(varargin,
                         "{?}",
                         mlfCreateColonIndex());
    }

  double evaluate_mx(mxArray* x_mx)
  {
    DOTRACE("evaluate_mx");

    mxArray* mx =  mlfFeval(mclValueVarargout(),
                            itsFunfcn,
                            x_mx,
                            getref(itsVarargin_ref),
                            NULL);
    double result = mxGetScalar(mx);
    mxDestroyArray(mx);
    return result;
  }

  virtual double doEvaluate(const Mtx& x)
  {
    return evaluate_mx(x.makeMxArray());
  }

public:
  MatlabFunction(mxArray* funfcn_mx, mxArray* varargin_mx) :
    MultivarFunction(),
    itsFunfcn(funfcn_mx),
    itsVarargin_ref(varargin_mx)
  {}

  virtual ~MatlabFunction() {}
};

int extractMaxIters(const mxArray* arr, int numModelParams)
{
  if (mxIsChar(arr))
    {
      if (MxWrapper::extractString(arr) == "200*numberofvariables")
        {
          return 200*numModelParams;
        }
      else
        {
          mexErrMsgTxt("Option must be an integer value "
                       "if not the default.");
        }
    }

  return int(mxGetScalar(arr));
}

namespace
{
  MexBuf* mexBuf = 0;

  std::streambuf* coutOrigBuf = 0;
  std::streambuf* cerrOrigBuf = 0;
}

void InitializeModule_doSimplex(void)
{
  mexPrintf("loading doSimplex mex file\n");

  mexBuf = new MexBuf;
#ifdef MIPS_PRO
  std::cout = mexBuf;
  std::cerr = mexBuf;
#else
  coutOrigBuf = std::cout.rdbuf(mexBuf);
  cerrOrigBuf = std::cerr.rdbuf(mexBuf);
#endif
}

void TerminateModule_doSimplex(void)
{
  mexPrintf("unloading doSimplex mex file\n");

  Util::Prof::printAtExit(false);

  // For some reason we get crashses if we try to reset these to their
  // original streambuf*'s; setting them to 0 effectively just shuts the
  // streams down
  std::cout.rdbuf(0);
  std::cerr.rdbuf(0);

  delete mexBuf;
}

static mxArray * MdoSimplex(mxArray * * fval,
                            mxArray * * exitflag,
                            mxArray * * output,
                            int nargout_,
                            mxArray * funfcn,
                            mxArray * x_in,
                            mxArray * printtype,
                            mxArray * tolx,
                            mxArray * tolf,
                            mxArray * maxfun,
                            mxArray * maxiter,
                            mxArray * debugFlags,
                            mxArray * varargin);

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
  mxArray * mprhs[9];
  mxArray * mplhs[4];
  int i;
  if (nlhs > 4)
    {
      mexErrMsgTxt("Run-time Error: File: doSimplex Line: 1 Column: 1 "
                   "The function \"doSimplex\" was called with more "
                   "than the declared number of outputs (4).");
    }
  for (i = 0; i < 4; ++i)
    {
      mplhs[i] = mclGetUninitializedArray();
    }
  for (i = 0; i < 8 && i < nrhs; ++i)
    {
      mprhs[i] = prhs[i];
    }
  for (; i < 8; ++i)
    {
      mprhs[i] = NULL;
    }
  mlfEnterNewContext(0,
                     8,
                     mprhs[0],
                     mprhs[1],
                     mprhs[2],
                     mprhs[3],
                     mprhs[4],
                     mprhs[5],
                     mprhs[6],
                     mprhs[7]);
  mprhs[8] = NULL;
  mlfAssign(&mprhs[8], mclCreateVararginCell(nrhs - 8, prhs + 8));

  mplhs[0] = MdoSimplex(&mplhs[1],
                        &mplhs[2],
                        &mplhs[3],
                        nlhs,
                        mprhs[0],
                        mprhs[1],
                        mprhs[2],
                        mprhs[3],
                        mprhs[4],
                        mprhs[5],
                        mprhs[6],
                        mprhs[7],
                        mprhs[8]);

  mlfRestorePreviousContext(0,
                            8,
                            mprhs[0],
                            mprhs[1],
                            mprhs[2],
                            mprhs[3],
                            mprhs[4],
			    mprhs[5],
                            mprhs[6],
                            mprhs[7]);

  plhs[0] = mplhs[0];
  for (i = 1; i < 4 && i < nlhs; ++i)
    {
      plhs[i] = mplhs[i];
    }
  for (; i < 4; ++i)
    {
      mxDestroyArray(mplhs[i]);
    }
  mxDestroyArray(mprhs[8]);
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
                            int nargout_,
                            mxArray * funfcn_mx,
                            mxArray * x_in,
                            mxArray * printtype_mx,
                            mxArray * tolx_mx,
                            mxArray * tolf_mx,
                            mxArray * maxfun_mx,
                            mxArray * maxiter_mx,
                            mxArray * debugFlags_mx,
                            mxArray * varargin)
{
DOTRACE("MdoSimplex");

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
  if (debugFlags_mx && mxIsStruct(debugFlags_mx))
    {
      mxArray* debugFlag = mxGetField(debugFlags_mx, 0, "debugFlag");

      if (debugFlag)
        {
          if (mxGetScalar(debugFlag) == -1)
            {
              Util::Prof::printAllProfData(std::cerr);
              return mxCreateScalarDouble(-1.0);
            }

          if (mxGetScalar(debugFlag) == -2)
            {
              Util::Prof::resetAllProfData();
              return mxCreateScalarDouble(-2.0);
            }
        }
    }
#endif

  try
    {
      // numModelParams = prod(size(x));
      const int numModelParams = mxGetM(x_in) * mxGetN(x_in);

      // Convert to inline function as needed.

      // FIXME Since this requires "object-oriented" programming, we can't
      // keep this and still use the MATLAB compiler

      // %funfcn = fcnchk(funfcn,length(varargin));

      MatlabFunction objective(funfcn_mx, varargin);

      SimplexOptimizer opt(objective,
                           Mtx(x_in),
                           MxWrapper::extractString(printtype_mx),
                           numModelParams,
                           extractMaxIters(maxfun_mx, numModelParams),
                           extractMaxIters(maxiter_mx, numModelParams),
                           mxGetScalar(tolx_mx),
                           mxGetScalar(tolf_mx)
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

static const char vcid_doSimplex_cc[] = "$Header$";
#endif // !DOSIMPLEX_CC_DEFINED
