///////////////////////////////////////////////////////////////////////
//
// annealHelper.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Mon Feb 18 16:30:29 2002
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

#include "matlabfunction.h"
#include "mexbuf.h"
#include "mtx.h"
#include "mxwrapper.h"
#include "rutil.h"
#include "simplexoptimizer.h"

#include "util/error.h"
#include "util/strings.h"

#include <iostream>
#include <libmatlbm.h>

#include "util/trace.h"
#include "util/debug.h"

namespace Mx
{
  inline mxArray* getField(mxArray* structArray, const char* fieldName)
  {
    return MxWrapper::extractStructField(structArray, fieldName);
  }
}

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

  Mtx matlabRand(int mrows, int ncols)
  {
    mxArray* arr = mlfRand(mxCreateScalarDouble(mrows),
                           mxCreateScalarDouble(ncols), NULL);

    Mtx result(arr, Mtx::COPY);
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

struct Astate
{
  Astate(mxArray* arr, int numruns) :
    talking(mxGetScalar(Mx::getField(arr, "talk")) != 0.0),
    canUseMatrix(mxGetScalar(Mx::getField(arr, "canUseMatrix")) != 0.0),
    doNewton(mxGetScalar(Mx::getField(arr, "newton")) != 0.0),
    numTemps(int(mxGetScalar(Mx::getField(arr, "numTemps")))),
    tempRepeats(Mx::getField(arr, "x"), Mtx::COPY),
    tempScales(Mx::getField(arr, "tempScales"), Mtx::COPY),
    numModelParams(int(mxGetScalar(Mx::getField(arr, "numModelParams")))),
    numStartingPoints(int(mxGetScalar(Mx::getField(arr, "numStartingPoints")))),
    valueScalingRange(Mx::getField(arr, "valueScalingRange"), Mtx::COPY),
    bounds(Mx::getField(arr, "bounds"), Mtx::COPY),
    deltas(Mx::getField(arr, "currentDeltas"), Mtx::COPY),
    numFunEvals(numruns, 1),
    bestCost(numruns, 1),
    modelHist(numModelParams, int(tempRepeats.sum())),
    mhat(numModelParams, numruns),
    energy(int(tempRepeats.sum()), numruns),
    startValues(Mx::getField(arr, "centers"), Mtx::COPY)
  {
    energy.setAll(std::numeric_limits<double>::max());
  }

  const bool talking;
  const bool canUseMatrix;
  const bool doNewton;
  const int numTemps;
  const Mtx tempRepeats;
  const Mtx tempScales;
  const int numModelParams;
  const int numStartingPoints;
  const Mtx valueScalingRange;
  const Mtx bounds;
  Mtx deltas;
  Mtx numFunEvals;
  Mtx bestCost;
  Mtx modelHist;
  Mtx mhat;
  Mtx energy;
  Mtx startValues;
};

//---------------------------------------------------------------------
//
// class AnnealingRun
//
//---------------------------------------------------------------------

class AnnealingRun
{
private:
  Astate& itsAstate;
  int itsRunNum;
  int itsNvisits;
  double itsCriticalTemp;
  Mtx itsMinUsedParams;
  Mtx itsMaxUsedParams;
  fstring itsFunFunName;
  int itsNvararg;
  mxArray** itsPvararg;

public:
  AnnealingRun(const fstring& funcName,
               Astate& astate,
               int nvararg,
               mxArray** pvararg)
    :
    itsAstate(astate),
    itsRunNum(0),
    itsNvisits(0),
    itsCriticalTemp(std::numeric_limits<double>::max()),
    itsMinUsedParams(0,0),
    itsMaxUsedParams(0,0),
    itsFunFunName(funcName),
    itsNvararg(nvararg),
    itsPvararg(pvararg)
  {}

  void oneRun();

  // Returns the cost at the best model that was found
  double visitParameters(Mtx& bestModel, const double temp);

  Mtx doFuncEvals(const Mtx& models)
  {
    if (itsAstate.canUseMatrix)
      return doParallelFuncEvals(models, itsFunFunName,
                                 itsNvararg, itsPvararg);

    return doSerialFuncEvals(models, itsFunFunName,
                             itsNvararg, itsPvararg);
  }

  Mtx makeRandomModels()
  {
    Mtx models = matlabRand(itsAstate.numModelParams,
                            itsAstate.numStartingPoints);

    models -= 0.5;
    models *= 2;

    for (int r = 0; r < models.mrows(); ++r)
      {
        models.row(r) *= itsAstate.deltas.at(r);
        models.row(r) += itsAstate.startValues.at(r);
      }

    return models;
  }

  Mtx createStartingModel()
  {
    DOTRACE("createStartingModel");

    const Mtx startingModels = makeRandomModels();

    const Mtx startingCosts = doFuncEvals(startingModels);

    itsCriticalTemp =
      log10(startingCosts.sum() / startingCosts.nelems())
      - itsAstate.tempScales.at(itsRunNum);

    int startingPos = 0;
    const double startingCost = startingCosts.min(&startingPos);

    printRunHeader(startingCost);

    Mtx startingModel = startingModels.column(startingPos);

    itsMinUsedParams = startingModel;
    itsMaxUsedParams = startingModel;

    return startingModel;
  }

  void printRunHeader(double startingCost)
  {
    DOTRACE("printRunHeader");

    if (!itsAstate.talking) return;

    mexPrintf("\nStarting cost %7.2f", startingCost);
    mexPrintf("\n\nBeginning run #%02d. Critical temperature at %3.2f.\n",
              itsRunNum+1, pow(10.0, itsCriticalTemp));
    mexPrintf("------------------------------------------------\n\n");
    mexPrintf("f-Calls\t\tTemperature\tMinimum f-Value\n");
    mexPrintf("------------------------------------------------\n");
  }

  void displayParams(const Mtx& model, double cost)
  {
    if (!itsAstate.talking) return;

    mexPrintf("\nparams: ");
    for (int i = 0; i < model.nelems(); ++i)
      mexPrintf("%f ", model.at(i));
    mexPrintf("\ncost: %7.4f", cost);
    mexPrintf("\n");
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
    for (int i = 0; i < itsAstate.deltas.nelems(); ++i)
      {
        itsAstate.deltas.at(i) =
          0.75 * (itsMaxUsedParams.at(i) - itsMinUsedParams.at(i));
      }
  }

  void updateBests()
  {
    // FIXME ought to use a smarter algorithm to keep track of the best cost

    Mtx currentEnergy = itsAstate.energy.column(itsRunNum);

    int best_pos = 0;
    const double best_energy = currentEnergy.min(&best_pos);

    itsAstate.bestCost.at(itsRunNum) = best_energy;
    itsAstate.mhat.column(itsRunNum) = itsAstate.modelHist.column(best_pos);

    displayParams(itsAstate.mhat.column(itsRunNum), itsAstate.bestCost.at(itsRunNum));
  }

  void runSimplex()
  {
    mxArray* funfun_mx = 0;
    mlfAssign(&funfun_mx, mxCreateString(itsFunFunName.c_str()));

    MatlabFunction objective(funfun_mx, itsNvararg, itsPvararg);

    SimplexOptimizer opt(objective,
                         Mtx(itsAstate.mhat.column(itsRunNum)),
                         fstring("notify"),
                         itsAstate.numModelParams,
                         10000000, // maxFunEvals
                         10000, // maxIter
                         1e-4, // tolx
                         1e-4 // tolf
                         );

    /*int exitFlag =*/ opt.optimize();

    double Ostar = opt.bestFval();

    Mtx mstar = opt.bestParams();

    if (Ostar < itsAstate.bestCost.at(itsRunNum))
      {
        displayParams(mstar, Ostar);

        if (itsAstate.talking)
          mexPrintf("%d iterations\n", opt.iterCount());

        bool inBounds = true;

        for (int i = 0; i < mstar.nelems(); ++i)
          {
            double val = mstar.at(i);
            if (val < itsAstate.bounds.at(i, 0)) { inBounds = false; break; }
            if (val > itsAstate.bounds.at(i, 1)) { inBounds = false; break; }
          }

        if (inBounds)
          {
            itsAstate.mhat.column(itsRunNum) = mstar;
            itsAstate.bestCost.at(itsRunNum) = Ostar;
            if (itsAstate.talking)
              mexPrintf("\nSimplex method lowered cost "
                        "and remained within constraints.\n\n");
          }
        else
          {
            if (itsAstate.talking)
              mexPrintf("\nSimplex method lowered cost "
                        "but failed to remain within constraints.\n\n");
          }
      }


    mxDestroyArray(funfun_mx);
  }
};

double AnnealingRun::visitParameters(Mtx& bestModel, const double temp)
{
  int nevals = 0;
  int s = -1;

  Mtx costs(0,0);

  // for x = find(deltas' ~= 0)
  for (int x = 0; x < itsAstate.deltas.nelems(); ++x)
    {
      const double delta = itsAstate.deltas.at(x);

      if (delta == 0.0) continue;

      Mtx modelmatrix = makeTestModels(x,
                                       bestModel,
                                       itsAstate.valueScalingRange,
                                       delta,
                                       itsAstate.bounds);

      costs = doFuncEvals(modelmatrix);

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

  itsAstate.numFunEvals.at(itsRunNum) += nevals;

  return costs.at(s);
}

//---------------------------------------------------------------------
//
// AnnealingRun::oneRun()
//
//---------------------------------------------------------------------

void AnnealingRun::oneRun()
{
DOTRACE("AnnealingRun::oneRun");

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

  itsNvisits = 0;
  itsCriticalTemp = std::numeric_limits<double>::max();
  itsMinUsedParams = Mtx(0,0);
  itsMaxUsedParams = Mtx(0,0);

  Mtx bestModel = createStartingModel();

  for (int temps_i = 0; temps_i < itsAstate.numTemps; ++temps_i)
    {
      // so the temperatures range from 10^(crit_temp+1) ... 10^(crit_temp-1)
      const double temp =
        pow(10.0, (itsCriticalTemp + 1.0
                   - 2.0*double(temps_i)/(itsAstate.numTemps-1)));

      const int temp_repeat = int(itsAstate.tempRepeats.at(temps_i));

      for (int repeat = 0; repeat < temp_repeat; ++repeat)
        {
          ++itsNvisits;

          if (itsAstate.talking && (itsNvisits % 10 == 0))
            {
              mexPrintf("%7d\t\t%7.2f\t\t%7.2f\n",
                        int(itsAstate.numFunEvals.at(itsRunNum)),
                        temp,
                        itsAstate.energy.column(itsRunNum).min());
            }

          itsAstate.energy.at(itsNvisits-1,itsRunNum) =
            visitParameters(bestModel, temp);

          updateUsedParams(bestModel);

          itsAstate.modelHist.column(itsNvisits-1) = bestModel;
        }
    }

  updateDeltas();

  updateBests();

  if (itsAstate.doNewton) runSimplex();

  itsAstate.startValues.column(0) = itsAstate.mhat.column(itsRunNum);

  ++itsRunNum;
}

///////////////////////////////////////////////////////////////////////
//
// class Annealer
//
///////////////////////////////////////////////////////////////////////

class Annealer
{
public:

  mxArray* go(mxArray* astate_mx,
              const fstring& func_name,
              int nvararg,
              mxArray** pvararg)
  {
    int numruns = int(mxGetScalar(Mx::getField(astate_mx, "numruns")));

    Astate astate(astate_mx, numruns);

    AnnealingRun ar(func_name, astate, nvararg, pvararg);

    for (int i = 0; i < numruns; ++i)
      {
        ar.oneRun();
      }

    const char* fieldnames[] =
      {
        "bestCost",
        "mhat",
        "model",
        "energy",
        "numFunEvals"
      };

    mxArray* output = mxCreateStructMatrix(1, 1, 5, fieldnames);

    mxSetField(output, 0, "bestCost", astate.bestCost.makeMxArray());
    mxSetField(output, 0, "mhat", astate.mhat.makeMxArray());
    mxSetField(output, 0, "model", astate.modelHist.makeMxArray());
    mxSetField(output, 0, "energy", astate.energy.makeMxArray());
    mxSetField(output, 0, "numFunEvals", astate.numFunEvals.makeMxArray());

    return output;
  }
};

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
      Annealer a;

      plhs[0] = a.go(prhs[0], // astate
                     MxWrapper::extractString(prhs[1]), // funcName
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

  mlfRestorePreviousContext(0, NDECLARED, prhs[0], prhs[1]);
}

static const char vcid_annealHelper_cc[] = "$Header$";
#endif // !ANNEALHELPER_CC_DEFINED
