///////////////////////////////////////////////////////////////////////
//
// annealVisitParameters.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Fri Feb 15 10:26:09 2002
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
// "annealVisitParameters"
//
///////////////////////////////////////////////////////////////////////

#ifndef ANNEALVISITPARAMETERS_CC_DEFINED
#define ANNEALVISITPARAMETERS_CC_DEFINED

#include "annealVisitParameters.h"

#include "mexbuf.h"
#include "mtx.h"
#include "mxwrapper.h"
#include "rutil.h"

#include "util/error.h"
#include "util/strings.h"

#include <iostream>
#include <libmatlbm.h>

#include "util/trace.h"
#include "util/debug.h"

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

  double matlabRand()
  {
    mxArray* arr = mlfNRand(1, NULL);
    double result = mxGetScalar(arr);
    mxDestroyArray(arr);
    return result;
  }
}

#ifdef MXDEBUG
#define mxDestroyArray(x)  debugDestroyArray(x, #x, __FILE__, __LINE__);
#endif

void InitializeModule_annealVisitParameters()
{
  mexPrintf("loading 'annealVisitParameters mex file\n");

  mexBuf = new MexBuf;
#ifdef MIPS_PRO
  std::cout = mexBuf;
  std::cerr = mexBuf;
#else
  coutOrigBuf = std::cout.rdbuf(mexBuf);
  cerrOrigBuf = std::cerr.rdbuf(mexBuf);
#endif
}

void TerminateModule_annealVisitParameters()
{
  mexPrintf("unloading 'annealVisitParameters' mex file\n");

  Util::Prof::printAtExit(false);

  // For some reason we get crashses if we try to reset these to their
  // original streambuf*'s; setting them to 0 effectively just shuts the
  // streams down
  std::cout.rdbuf(0);
  std::cerr.rdbuf(0);

  delete mexBuf;
}

//---------------------------------------------------------------------
//
// makeTestModels()
//
//---------------------------------------------------------------------

Mtx makeTestModels(int x,
                   const Mtx& bestModel,
                   const Mtx& valueScalingRange,
                   const double delta,
                   const Mtx& bounds)
{
DOTRACE("makeTestModels");

  const double current_x_val = bestModel.at(x);

  const double lowerbound_x = bounds.at(x,0);
  const double upperbound_x = bounds.at(x,1);

  int num_testmodels = 0;
  {for (MtxConstIter iter = valueScalingRange.rowIter(0); iter.hasMore(); ++iter)
    {
      double value = current_x_val + (*iter) * delta;
      if ( (value > lowerbound_x) && (value < upperbound_x) )
        ++num_testmodels;
    }
  }

  Mtx newModels(bestModel.mrows(), num_testmodels);

  const Slice bestModelCol = bestModel.column(0);

  for (int i = 0; i < num_testmodels; ++i)
    {
      newModels.column(i) = bestModelCol;
    }

  {
    int col = 0;
    for (MtxConstIter iter = valueScalingRange.rowIter(0); iter.hasMore(); ++iter)
    {
      double value = current_x_val + (*iter) * delta;
      if ( (value > lowerbound_x) && (value < upperbound_x) )
        newModels.at(x, col++) = value;
    }
  }

  return newModels;
}

//---------------------------------------------------------------------
//
// doFuncEvals()
//
//---------------------------------------------------------------------

Mtx doParallelFuncEvals(const Mtx& models,
                        mxArray* func,
                        mxArray* /*varargin_mx*/,
                        int nvararg,
                        mxArray** pvararg)
{
DOTRACE("doParallelFuncEvals");

  mxArray* costs_mx = 0;
  mxArray* models_mx = 0; mlfAssign(&models_mx, models.makeMxArray());
  mclCopyArray(&func);
//   mclCopyArray(&varargin_mx);

  // costs = feval(func, models, varargin{:});
#if 0
  mlfAssign(&costs_mx,
            mlfFeval(mclValueVarargout(),
                     func,
                     models_mx,
                     mlfIndexRef(varargin_mx, "{?}", mlfCreateColonIndex()),
                     NULL));
#else

  mxArray* plhs[1] = { 0 };

  const int MAX_NRHS = 32;
  mxArray* prhs[MAX_NRHS];

  prhs[0] = models_mx;
  // need this? seems like no. prhs[0] = mxDuplicateArray(models_mx);

  int nrhs = 1;

//   const int nvararg = mxGetNumberOfElements(varargin_mx);

  for (int i = 0; i < nvararg && nrhs < MAX_NRHS; ++i)
    {
      DebugEvalNL(i);
//       prhs[nrhs++] = mxDuplicateArray(mxGetCell(varargin_mx, i));
      prhs[nrhs++] = mxDuplicateArray(pvararg[i]);
    }

  fstring cmd_name = MxWrapper::extractString(func);

  int result = mexCallMATLAB(1, plhs, nrhs, prhs, cmd_name.c_str());

  if (result != 0) mexErrMsgTxt("mexCallMATLAB failed in doFuncEvals");

  mlfAssign(&costs_mx, plhs[0]);
#endif

  mclValidateOutput(costs_mx, 1, 1, "costs", "annealVisitParameters/doFuncEvals");

//   mxDestroyArray(varargin_mx);
  mxDestroyArray(func);
  mxDestroyArray(models_mx);

  Mtx costs(costs_mx, Mtx::COPY);

  mxDestroyArray(costs_mx);

  return costs;
}

Mtx doSerialFuncEvals(const Mtx& models,
                      mxArray* func,
                      mxArray* varargin_mx)
{
DOTRACE("doSerialFuncEvals");

  mxArray* models_mx = 0; mlfAssign(&models_mx, models.makeMxArray());
  mclCopyArray(&func);
  mclCopyArray(&varargin_mx);

  const int NM = models.ncols();

  Mtx costs(NM, 1);

  for (int e = 0; e < NM; ++e)
    {
      //costs(e) = feval(func, models(:,e), varargin{:});
      mxArray* result = mlfFeval(mclValueVarargout(),
                                 func,
                                 mclArrayRef2(models_mx,
                                              mlfCreateColonIndex(),
                                              mlfScalar(e+1)),
                                 mlfIndexRef(varargin_mx,
                                             "{?}",
                                             mlfCreateColonIndex()),
                                 NULL);

      costs.at(e) = mxGetScalar(result);

      mxDestroyArray(result);
    }

  mxDestroyArray(varargin_mx);
  mxDestroyArray(func);
  mxDestroyArray(models_mx);

  return costs;
}

Mtx doFuncEvals(bool canUseMatrix,
                const Mtx& models,
                mxArray* func,
                mxArray* varargin_mx,
                int nvararg,
                mxArray** pvararg)
{
  if (canUseMatrix)
    {
      return doParallelFuncEvals(models, func, varargin_mx, nvararg, pvararg);
    }

  return doSerialFuncEvals(models, func, varargin_mx);
}

//---------------------------------------------------------------------
//
// makePDF()
//
//---------------------------------------------------------------------

Mtx makePDF(const Mtx& temp, const Mtx& costs)
{
DOTRACE("makePDF");

  for (int i = 0; i < costs.nelems(); ++i)
    if (isnan(costs.at(i)))
      {
        mexErrMsgTxt("Warning: the cost function generated a NaN "
                     "for one or more models.");
      }

  // Scales cost vector and calculates exponential probability distribution.
  // The scaling is necessary to permit wide ranges in temperature.  Internal
  // function for simulated annealing algorithm.

  const double toobig = 708.3964185322641;

  Mtx pdf = costs; pdf /= temp.at(0);

  const double mpdf = pdf.max();

  if (mpdf > toobig)
    {
      const double scale = mpdf/toobig;

      pdf *= (-1.0/scale); pdf.apply(exp);

      pdf /= pdf.max();

      pdf.applyF(ToPow(scale));
    }
  else
    {
      pdf *= -1; pdf.apply(exp);

      pdf /= pdf.max();
    }

  pdf /= pdf.sum();

  return pdf;
}

//---------------------------------------------------------------------
//
// sampleFromPdf()
//
//---------------------------------------------------------------------

int sampleFromPdf(const Mtx& temp, const Mtx& costs)
{
DOTRACE("sampleFromPdf");

  Mtx dist = makePDF(temp, costs);

  const double cutoff = matlabRand();

  double cumsum = 0.0;
  int s = -1;
  for (int i = 0; i < dist.nelems(); ++i)
    {
      cumsum += dist.at(i);

      if (cumsum >= cutoff) { s = i; break; }
    }

  if (s < 0)
    {
      mexErrMsgTxt("snafu in sampleFromPdf");
    }

  return s;
}

/*
 * The function "MannealVisitParameters" is the implementation version of the
 * "annealVisitParameters" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 1-28). It contains the actual compiled code for that M-function. It
 * is a static function and must only be called from one of the interface
 * functions, appearing below.
 */

mxArray* MannealVisitParameters(int nargout_,
                                mxArray* bestModel_mx,
                                mxArray* valueScalingRange_mx,
                                mxArray* deltas_mx,
                                mxArray* bounds_mx,
                                mxArray* canUseMatrix_mx,
                                mxArray* FUN_mx,
                                mxArray* temp_mx,
                                mxArray* varargin_mx,
                                int nvararg,
                                mxArray** pvararg)
{
DOTRACE("MannealVisitParameters");

  try
    {

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
      if (varargin_mx && mxGetScalar(mxGetCell(varargin_mx,0)) == -1)
        {
          Util::Prof::printAllProfData(std::cerr);
          return mxCreateScalarDouble(-1.0);
        }

      if (varargin_mx && mxGetScalar(mxGetCell(varargin_mx,0)) == -2)
        {
          Util::Prof::resetAllProfData();
          return mxCreateScalarDouble(-2.0);
        }
#endif

      mclCopyArray(&bestModel_mx);

      Mtx bestModel(bestModel_mx, Mtx::REFER);

      const bool canUseMatrix = (mxGetPr(canUseMatrix_mx)[0] != 0.0);

      int nevals = 0;
      int s_zerobased = -1;

      const Mtx valueScalingRange(valueScalingRange_mx, Mtx::BORROW);
      const Mtx deltas(deltas_mx, Mtx::BORROW);
      const Mtx bounds(bounds_mx, Mtx::BORROW);

      Mtx costs(0,0);

      // for x = find(deltas' ~= 0)
      for (int x = 0; x < deltas.nelems(); ++x)
        {
          const double delta = deltas.at(x);

          if (delta == 0.0) continue;

          Mtx modelmatrix = makeTestModels(x,
                                           bestModel,
                                           valueScalingRange,
                                           delta,
                                           bounds);

          // costs = doFuncEvals(canUseMatrix, modelmatrix, FUN, varargin{:});
          costs = doFuncEvals(canUseMatrix,
                              modelmatrix,
                              FUN_mx,
                              mlfIndexRef(varargin_mx,
                                          "{?}",
                                          mlfCreateColonIndex()),
                              nvararg,
                              pvararg);

          // S.nevals = S.nevals + length(costs);
          nevals += costs.nelems();;

          // Sample from probability distribution
          s_zerobased = sampleFromPdf(Mtx(temp_mx, Mtx::BORROW), costs);

          bestModel.at(x, 0) = modelmatrix.at(x, s_zerobased);
        }

      if (s_zerobased < 0)
        {
          mexErrMsgTxt("didn't run any annealing iterations because "
                       "there were no non-zero deltas");
        }

      const char* fieldNames[] = { "nevals", "newModel", "cost" };
      mxArray* output = mxCreateStructMatrix(1,1,3,fieldNames);

      // S.nevals = 0;
      mxSetField(output, 0, "nevals", mxCreateScalarDouble(nevals));

      // S.newModel = bestModel;
      mxSetField(output, 0, "newModel", bestModel_mx);

      // S.cost = costs(s);
      mxSetField(output, 0, "cost",
                 mxCreateScalarDouble(costs.at(s_zerobased)));

      return output;
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
 * The function "mlxAnnealVisitParameters" contains the feval interface for the
 * "annealVisitParameters" M-function from file
 * "/cit/rjpeters/science/psyphy/classmodels/matlab/annealVisitParameters.m"
 * (lines 1-28). The feval function calls the implementation version of
 * annealVisitParameters through this function. This function processes any
 * input arguments and passes them to the implementation version of the
 * function, appearing above.
 */
void mlxAnnealVisitParameters(int nlhs, mxArray* plhs[],
                              int nrhs, mxArray* prhs[])
{
DOTRACE("mlxAnnealVisitParameters");

  if (nlhs > 1)
    {
      mexErrMsgTxt("Error: annealVisitParameters was called with more "
                   "than the declared number of outputs (1).");
    }
  if (nrhs < 7)
    {
      mexErrMsgTxt("Error: annealVisitParameters was called with fewer "
                   "than the declared number of inputs (7).");
    }

  mxArray* varargin = NULL;

  mlfEnterNewContext(0,
                     7,
                     prhs[0],
                     prhs[1],
                     prhs[2],
                     prhs[3],
                     prhs[4],
                     prhs[5],
                     prhs[6]);
  varargin = NULL;
  mlfAssign(&varargin, mclCreateVararginCell(nrhs - 7, prhs + 7));

  int nvararg = nrhs - 7;
  mxArray** pvararg = prhs + 7;

  plhs[0] = MannealVisitParameters(nlhs,
                                   prhs[0],
                                   prhs[1],
                                   prhs[2],
                                   prhs[3],
                                   prhs[4],
                                   prhs[5],
                                   prhs[6],
                                   varargin,
                                   nvararg,
                                   pvararg);

  mlfRestorePreviousContext(0,
                            7,
                            prhs[0],
                            prhs[1],
                            prhs[2],
                            prhs[3],
                            prhs[4],
                            prhs[5],
                            prhs[6]);

  mxDestroyArray(varargin);
}

static const char vcid_annealVisitParameters_cc[] = "$Header$";
#endif // !ANNEALVISITPARAMETERS_CC_DEFINED
