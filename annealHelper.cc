///////////////////////////////////////////////////////////////////////
//
// annealVisitParameters.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Fri Feb 15 15:35:31 2002
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
// doParallelFuncEvals()
//
//---------------------------------------------------------------------

Mtx doParallelFuncEvals(const Mtx& models,
                        const fstring& func_name,
                        int nvararg,
                        mxArray** pvararg)
{
DOTRACE("doParallelFuncEvals");

  mxArray* costs_mx = 0;
  mxArray* models_mx = 0; mlfAssign(&models_mx, models.makeMxArray());

  const int MAX_NRHS = 32;
  mxArray* prhs[MAX_NRHS];

  // We don't need to call mxDuplicateArray, since models_mx is a bound
  // variable due to the mlfAssign() above; that's also why we have to
  // explicitly destroy it later on
  prhs[0] = models_mx;

  int nrhs = 1;

  for (int i = 0; i < nvararg && nrhs < MAX_NRHS; ++i)
    {
      prhs[nrhs++] = mxDuplicateArray(pvararg[i]);
    }

  // costs = feval(func, models, varargin{:});
  int result = mexCallMATLAB(1, &costs_mx, nrhs, prhs, func_name.c_str());

  if (result != 0) mexErrMsgTxt("mexCallMATLAB failed in doFuncEvals");

  Mtx costs(costs_mx, Mtx::COPY);

  mxDestroyArray(models_mx);
  mxDestroyArray(costs_mx);

  return costs;
}

//---------------------------------------------------------------------
//
// doSerialFuncEvals()
//
//---------------------------------------------------------------------

Mtx doSerialFuncEvals(const Mtx& models,
                      const fstring& func_name,
                      int nvararg,
                      mxArray** pvararg)
{
DOTRACE("doSerialFuncEvals");

  mxArray* models_mx = 0; mlfAssign(&models_mx, models.makeMxArray());

  const int NM = models.ncols();

  Mtx costs(NM, 1);

  mxArray* plhs[1] = { 0 };

  const int MAX_NRHS = 32;
  mxArray* prhs[MAX_NRHS];

  prhs[0] = 0;

  int nrhs = 1;

  for (int i = 0; i < nvararg && nrhs < MAX_NRHS; ++i)
    {
      prhs[nrhs] = 0; mlfAssign(prhs+nrhs, mxDuplicateArray(pvararg[i]));
      ++nrhs;
    }

  for (int e = 0; e < NM; ++e)
    {
      //costs(e) = feval(func, models(:,e), varargin{:});
      prhs[0] = mclArrayRef2(models_mx,
                             mlfCreateColonIndex(),
                             mlfScalar(e+1));

      int result = mexCallMATLAB(1, plhs, nrhs, prhs, func_name.c_str());

      if (result != 0) mexErrMsgTxt("mexCallMATLAB failed in doFuncEvals");

      costs.at(e) = mxGetScalar(plhs[0]);

      mxDestroyArray(plhs[0]);
    }

  while (--nrhs >= 1)
    {
      mxDestroyArray(prhs[nrhs]);
    }

  mxDestroyArray(models_mx);

  return costs;
}

//---------------------------------------------------------------------
//
// makePDF()
//
//---------------------------------------------------------------------

Mtx makePDF(double temp, const Mtx& costs)
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

  Mtx pdf = costs; pdf /= temp;

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

int sampleFromPdf(double temp, const Mtx& costs)
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

//---------------------------------------------------------------------
//
// annealVisitParameters()
//
//---------------------------------------------------------------------

namespace
{
  struct VisitResult
  {
    VisitResult(int n, mxArray* b, double c) :
      nevals(n), newModel(b), cost(c) {}

    int const nevals;
    mxArray* const newModel;
    double const cost;
  };
}

VisitResult annealVisitParameters(mxArray* bestModel_mx,
                                  const Mtx& valueScalingRange,
                                  const Mtx& deltas,
                                  const Mtx& bounds,
                                  const bool canUseMatrix,
                                  const fstring& func_name,
                                  const double temp,
                                  int nvararg,
                                  mxArray** pvararg)
{

  mclCopyArray(&bestModel_mx);

  Mtx bestModel(bestModel_mx, Mtx::REFER);

  int nevals = 0;
  int s = -1;

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
      if (canUseMatrix)
        costs = doParallelFuncEvals(modelmatrix, func_name,
                                    nvararg, pvararg);
      else
        costs = doSerialFuncEvals(modelmatrix, func_name,
                                  nvararg, pvararg);

      // S.nevals = S.nevals + length(costs);
      nevals += costs.nelems();

      // Sample from probability distribution
      s = sampleFromPdf(temp, costs);

      bestModel.at(x, 0) = modelmatrix.at(x, s);
    }

  if (s < 0)
    {
      mexErrMsgTxt("didn't run any annealing iterations because "
                   "there were no non-zero deltas");
    }

  return VisitResult(nevals, bestModel_mx, costs.at(s));
}

//---------------------------------------------------------------------
//
// annealHelper()
//
//---------------------------------------------------------------------

mxArray* annealHelper(mxArray* old_astate_mx,
                      const Mtx& valueScalingRange,
                      const Mtx& deltas,
                      const Mtx& bounds,
                      const bool canUseMatrix,
                      const fstring& func_name,
                      int nvararg,
                      mxArray** pvararg)
{
DOTRACE("annealHelper");

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
 if (nvararg > 0 && int(mxGetScalar(pvararg[0])) == -1)
   {
     Util::Prof::printAllProfData(std::cerr);
     return mxCreateScalarDouble(-1.0);
   }
 if (nvararg > 0 && int(mxGetScalar(pvararg[0])) == -2)
   {
     Util::Prof::resetAllProfData();
     return mxCreateScalarDouble(-2.0);
   }
#endif

  mxArray* const astate_mx = mxDuplicateArray(old_astate_mx);

  const int astate_t = int(mxGetScalar(mxGetField(astate_mx, 0, "t")));

  for (int w_onebased = 1; w_onebased <= astate_t; ++w_onebased)
    {
      const Mtx astate_T(mxGetField(astate_mx, 0, "T"), Mtx::BORROW);

      const double temp = astate_T.at(w_onebased-1);

      int astate_c = int(mxGetScalar(mxGetField(astate_mx, 0, "c")));
      ++astate_c;
      mxSetField(astate_mx, 0, "c", mxCreateScalarDouble(astate_c));

      const bool astate_talk =
        (mxGetScalar(mxGetField(astate_mx, 0, "talk")) != 0.0);

      Mtx numFunEvals (mxGetField(astate_mx, 0, "numFunEvals"), Mtx::REFER);
      Mtx energy      (mxGetField(astate_mx, 0, "energy"), Mtx::REFER);

      const int k_onebased = int(mxGetScalar(mxGetField(astate_mx, 0, "k")));

      if (astate_talk && (astate_c % 10 == 0))
        {
          mexPrintf("%7d\t\t%7.2f\t\t%7.2f\n",
                    int(numFunEvals.at(k_onebased -1)),
                    temp,
                    energy.column(k_onebased-1).leftmost(astate_c-1).min());
        }

      VisitResult vresult =
        annealVisitParameters(mxGetField(astate_mx, 0, "bestModel"),
                              valueScalingRange,
                              deltas,
                              bounds,
                              canUseMatrix,
                              func_name,
                              temp,
                              nvararg,
                              pvararg);

      mxSetField(astate_mx, 0, "bestModel", vresult.newModel);

      numFunEvals.at(k_onebased-1) += vresult.nevals;

      energy.at(astate_c-1,k_onebased-1) = vresult.cost;

      const Mtx bestModel(mxGetField(astate_mx, 0, "bestModel"), Mtx::REFER);

      // Update the used parameter ranges
      {
        Mtx minUsedParams(mxGetField(astate_mx, 0, "minUsedParams"), Mtx::REFER);
        Mtx maxUsedParams(mxGetField(astate_mx, 0, "maxUsedParams"), Mtx::REFER);

        if (minUsedParams.nelems() != bestModel.nelems()
            || maxUsedParams.nelems() != bestModel.nelems())
          {
            mexErrMsgTxt("size mismatch in annealHelper()");
          }

        for (int i = 0; i < bestModel.nelems(); ++i)
          {
            minUsedParams.at(i) = std::min(double(minUsedParams.at(i)),
                                           bestModel.at(i));

            maxUsedParams.at(i) = std::max(double(maxUsedParams.at(i)),
                                           bestModel.at(i));
          }
      }

      // Update the model history matrix
      {
        Mtx modelHist(mxGetField(astate_mx, 0, "model"), Mtx::REFER);

        modelHist.column(astate_c-1) = bestModel;
      }
    }

  return astate_mx;
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

  const int NDECLARED = 6;

  if (nrhs < NDECLARED)
    {
      mexErrMsgTxt("Error: annealVisitParameters was called with fewer "
                   "than the declared number of inputs (7).");
    }

  mlfEnterNewContext(0, NDECLARED,
                     prhs[0], prhs[1], prhs[2],
                     prhs[3], prhs[4], prhs[5]);

  int nvararg = nrhs - NDECLARED;
  mxArray** pvararg = prhs + NDECLARED;

  try
    {
      plhs[0] = annealHelper
        (prhs[0], // astate
         Mtx(prhs[1], Mtx::BORROW), // valueScalingRange
         Mtx(prhs[2], Mtx::BORROW), // deltas
         Mtx(prhs[3], Mtx::BORROW), // bounds
         (mxGetScalar(prhs[4]) != 0.0), // canUseMatrix
         MxWrapper::extractString(prhs[5]), // func_name
         nvararg,
         pvararg);
    }
  catch (Util::Error& err)
    {
      mexErrMsgTxt(err.msg_cstr());
    }
  catch (...)
    {
      mexErrMsgTxt("an unknown C++ exception occurred.");
    }

  mlfRestorePreviousContext(0, NDECLARED,
                            prhs[0], prhs[1], prhs[2],
                            prhs[3], prhs[4], prhs[5]);
}

static const char vcid_annealVisitParameters_cc[] = "$Header$";
#endif // !ANNEALVISITPARAMETERS_CC_DEFINED
