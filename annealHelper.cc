///////////////////////////////////////////////////////////////////////
//
// annealHelper.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Sun Feb 17 08:25:24 2002
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

void InitializeModule_annealHelper()
{
  mexPrintf("loading 'annealHelper mex file\n");

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

//---------------------------------------------------------------------
//
// class AnnealingRun
//
//---------------------------------------------------------------------

class AnnealingRun
{
public:
  AnnealingRun(mxArray* old_astate_mx,
               const fstring& funcName,
               int nvararg,
               mxArray** pvararg)
    :
    itsAstate_mx(mxDuplicateArray(old_astate_mx)),
    itsRunNum(int(mxGetScalar(mxGetField(itsAstate_mx, 0, "k")))-1),
    itsTalking(mxGetScalar(mxGetField(itsAstate_mx, 0, "talk")) != 0.0),
    itsNvisits(0),
    itsCriticalTemp(mxGetScalar(mxGetField(itsAstate_mx, 0, "crit_temp"))),
    itsNumTemps(int(mxGetScalar(mxGetField(itsAstate_mx, 0, "numTemps")))),
    itsTempRepeats(mxGetField(itsAstate_mx, 0, "x"), Mtx::BORROW),
    itsNumFunEvals(mxGetField(itsAstate_mx, 0, "numFunEvals"), Mtx::REFER),
    itsEnergy(mxGetField(itsAstate_mx, 0, "energy"), Mtx::REFER),
    itsMinUsedParams(mxGetField(itsAstate_mx, 0, "bestModel"), Mtx::COPY),
    itsMaxUsedParams(mxGetField(itsAstate_mx, 0, "bestModel"), Mtx::COPY),
    itsDeltas(mxGetField(itsAstate_mx, 0, "currentDeltas"), Mtx::REFER),
    itsValueScalingRange(mxGetField(itsAstate_mx, 0, "valueScalingRange"),
                         Mtx::BORROW),
    itsBounds(mxGetField(itsAstate_mx, 0, "bounds"), Mtx::BORROW),
    itsCanUseMatrix(mxGetScalar(mxGetField(itsAstate_mx, 0, "canUseMatrix"))
                    != 0.0),
    itsFunFunName(funcName),
    itsNvararg(nvararg),
    itsPvararg(pvararg)
  {}

  mxArray* go();

  VisitResult visitParameters(mxArray* bestModel_mx,
                              const double temp);

  void printRunHeader()
  {
    if (!itsTalking) return;

    const double startingCost =
      mxGetScalar(mxGetField(itsAstate_mx, 0, "startingCost"));

    mexPrintf("\nStarting cost %7.2f", startingCost);
    mexPrintf("\n\nBeginning run #%02d. Critical temperature at %3.2f.\n",
              itsRunNum+1, pow(10.0, itsCriticalTemp));
    mexPrintf("------------------------------------------------\n\n");
    mexPrintf("f-Calls\t\tTemperature\tMinimum f-Value\n");
    mexPrintf("------------------------------------------------\n");
  }

  void updateUsedParams(const Mtx& model)
  {
    for (int i = 0; i < model.nelems(); ++i)
      {
        itsMinUsedParams.at(i) = std::min(double(itsMinUsedParams.at(i)),
                                          model.at(i));

        itsMaxUsedParams.at(i) = std::max(double(itsMaxUsedParams.at(i)),
                                          model.at(i));
      }
  }

  void updateDeltas()
  {
    for (int i = 0; i < itsDeltas.nelems(); ++i)
      {
        itsDeltas.at(i) =
          0.75 * (itsMaxUsedParams.at(i) - itsMinUsedParams.at(i));
      }
  }

  void updateModelHistory(const Mtx& bestModel)
  {
    Mtx modelHist(mxGetField(itsAstate_mx, 0, "model"), Mtx::REFER);

    modelHist.column(itsNvisits-1) = bestModel;
  }

  void updateBests()
  {
    // FIXME ought to use a smarter algorithm to keep track of the best cost

    int best_pos = 0;
    double best_energy = itsEnergy.at(0, itsRunNum);
    for (int i = 1; i < itsEnergy.mrows(); ++i)
      {
        if (itsEnergy.at(i, itsRunNum) < best_energy)
          {
            best_pos = i;
            best_energy = itsEnergy.at(i, itsRunNum);
          }
      }

    Mtx bestCost(mxGetField(itsAstate_mx, 0, "bestCost"), Mtx::REFER);
    Mtx mhat(mxGetField(itsAstate_mx, 0, "mhat"), Mtx::REFER);

    const Mtx modelHist(mxGetField(itsAstate_mx, 0, "model"), Mtx::BORROW);

    bestCost.at(itsRunNum) = best_energy;
    mhat.column(itsRunNum) = modelHist.column(best_pos);
  }

private:
  mxArray* const itsAstate_mx;
  const int itsRunNum;
  const bool itsTalking;
  int itsNvisits;
  const double itsCriticalTemp;
  const int itsNumTemps;
  Mtx itsTempRepeats;
  Mtx itsNumFunEvals;
  Mtx itsEnergy;
  Mtx itsMinUsedParams;
  Mtx itsMaxUsedParams;
  Mtx itsDeltas;
  const Mtx itsValueScalingRange;
  const Mtx itsBounds;
  const bool itsCanUseMatrix;
  fstring itsFunFunName;
  int itsNvararg;
  mxArray** itsPvararg;
};

VisitResult AnnealingRun::visitParameters(mxArray* bestModel_mx,
                                          const double temp)
{

  mclCopyArray(&bestModel_mx);

  Mtx bestModel(bestModel_mx, Mtx::REFER);

  int nevals = 0;
  int s = -1;

  Mtx costs(0,0);

  // for x = find(deltas' ~= 0)
  for (int x = 0; x < itsDeltas.nelems(); ++x)
    {
      const double delta = itsDeltas.at(x);

      if (delta == 0.0) continue;

      Mtx modelmatrix = makeTestModels(x,
                                       bestModel,
                                       itsValueScalingRange,
                                       delta,
                                       itsBounds);

      // costs = doFuncEvals(canUseMatrix, modelmatrix, FUN, varargin{:});
      if (itsCanUseMatrix)
        costs = doParallelFuncEvals(modelmatrix, itsFunFunName,
                                    itsNvararg, itsPvararg);
      else
        costs = doSerialFuncEvals(modelmatrix, itsFunFunName,
                                  itsNvararg, itsPvararg);

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
// AnnealingRun::go()
//
//---------------------------------------------------------------------

mxArray* AnnealingRun::go()
{
DOTRACE("AnnealingRun::go");

#if defined(LOCAL_DEBUG) || defined(LOCAL_PROF)
  if (itsNvararg > 0 && int(mxGetScalar(itsPvararg[0])) == -1)
    {
      Util::Prof::printAllProfData(std::cerr);
      return mxCreateScalarDouble(-1.0);
    }
  if (itsNvararg > 0 && int(mxGetScalar(itsPvararg[0])) == -2)
    {
      Util::Prof::resetAllProfData();
      return mxCreateScalarDouble(-2.0);
    }
#endif

  printRunHeader();

  for (int temps_i = 0; temps_i < itsNumTemps; ++temps_i)
    {
      // so the temperatures range from 10^(crit_temp+1) ... 10^(crit_temp-1)
      const double temp =
        pow(10.0, itsCriticalTemp + 1.0 - 2.0*double(temps_i)/(itsNumTemps-1));

      const int temp_repeat = int(itsTempRepeats.at(temps_i));

      for (int repeat = 0; repeat < temp_repeat; ++repeat)
        {
          ++itsNvisits;

          if (itsTalking && (itsNvisits % 10 == 0))
            {
              mexPrintf("%7d\t\t%7.2f\t\t%7.2f\n",
                        int(itsNumFunEvals.at(itsRunNum)),
                        temp,
                        itsEnergy.column(itsRunNum).leftmost(itsNvisits-1).min());
            }

          VisitResult vresult =
            visitParameters(mxGetField(itsAstate_mx, 0, "bestModel"),
                            temp);

          mxSetField(itsAstate_mx, 0, "bestModel", vresult.newModel);

          itsNumFunEvals.at(itsRunNum) += vresult.nevals;

          itsEnergy.at(itsNvisits-1,itsRunNum) = vresult.cost;

          const Mtx bestModel(vresult.newModel, Mtx::REFER);

          updateUsedParams(bestModel);

          updateModelHistory(bestModel);
        }
    }

  updateDeltas();

  updateBests();

  return itsAstate_mx;
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
      AnnealingRun ar
        (prhs[0], // astate
         MxWrapper::extractString(prhs[1]), // funcName
         nvararg,
         pvararg);

      plhs[0] = ar.go();
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
