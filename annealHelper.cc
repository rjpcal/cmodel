///////////////////////////////////////////////////////////////////////
//
// annealHelper.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar 23 17:17:00 2001
// written: Tue Feb 19 09:53:50 2002
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
#include "multivarfunction.h"
#include "mxwrapper.h"
#include "rutil.h"
#include "simplexoptimizer.h"

#include "util/error.h"
#include "util/strings.h"

#include <iostream>
#include <iomanip>
#include <libmatlbm.h>

#include "util/trace.h"
#include "util/debug.h"

namespace Mx
{
  inline mxArray* getField(mxArray* structArray, const char* fieldName)
  {
    return MxWrapper::extractStructField(structArray, fieldName);
  }

  inline int getIntField(mxArray* structArray, const char* fieldName)
  {
    return int(mxGetScalar(Mx::getField(structArray, fieldName)));
  }

  inline bool getBoolField(mxArray* structArray, const char* fieldName)
  {
    return (mxGetScalar(Mx::getField(structArray, fieldName)) != 0.0);
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
        throw Util::Error("Warning: the cost function generated a NaN "
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
      throw Util::Error("snafu in sampleFromPdf");
    }

  return s;
}

///////////////////////////////////////////////////////////////////////
//
// AnnealOpts
//
///////////////////////////////////////////////////////////////////////

struct AnnealOpts
{
private:
  Mtx makeScalingRange(int gridSpacing)
  {
    Mtx result(1, 2*gridSpacing+1);

    for (int i = 0; i < gridSpacing; ++i)
      {
        double val = pow(2.0, -(i+1));
        result.at(i) = -val;
        result.at(2*gridSpacing-i) = val;
      }

    result.at(gridSpacing) = 0.0;

    return result;
  }

  Mtx makeCoolingSchedule(int scale)
  {
    const int N = 9;

    double repeats[N] = { 1.0, 2.0, 4.0, 6.0, 10.0, 6.0, 4.0, 2.0, 1.0 };

    Mtx result(1, N);

    for (int i = 0; i < N; ++i)
      {
        result.at(i) = scale * repeats[i];
      }

    return result;
  }

public:
  AnnealOpts(mxArray* arr) :
    bounds(Mx::getField(arr, "bounds"), Mtx::COPY),
    deltas(Mx::getField(arr, "deltas"), Mtx::COPY),
    centers(Mx::getField(arr, "centers"), Mtx::COPY),
    numRuns(Mx::getIntField(arr, "numruns")),
    talking(Mx::getBoolField(arr, "talk")),
    canUseMatrix(Mx::getBoolField(arr, "canUseMatrix")),
    doNewton(Mx::getBoolField(arr, "newton")),
    coolingSchedule(makeCoolingSchedule(Mx::getIntField(arr, "scale"))),
    numTemps(coolingSchedule.nelems()),
    tempScales(Mx::getField(arr, "tempScales"), Mtx::COPY),
    numModelParams(bounds.mrows()),
    numStartingPoints(std::max(100, 20*numModelParams)),
    valueScalingRange(makeScalingRange(Mx::getIntField(arr, "gridSpacing")))
  {
    if (bounds.ncols() != 2)
      {
        throw Util::Error("'bounds' must have two columns (lo and hi bound)");
      }

    for (int r = 0; r < bounds.mrows(); ++r)
      {
        if (bounds.at(r,0) > bounds.at(r,1))
          throw Util::Error("invalid 'bounds' (lo value was greater than hi)");
      }
  }

  const Mtx bounds;
  const Mtx deltas;
  const Mtx centers;
  const int numRuns;
  const bool talking;
  const bool canUseMatrix;
  const bool doNewton;
  const Mtx coolingSchedule;
  const int numTemps;
  const Mtx tempScales;
  const int numModelParams;
  const int numStartingPoints;
  const Mtx valueScalingRange;
};

//---------------------------------------------------------------------
//
// class AnnealingOptimizer
//
//---------------------------------------------------------------------

class AnnealingOptimizer
{
private:
  AnnealOpts& itsOpts;
  int itsRunNum;
  Mtx itsDeltas;
  Mtx itsNumFunEvals;
  Mtx itsBestCosts;
  Mtx itsModelHist; // FIXME does this really need to be a member variable?
  Mtx itsBestModels;
  Mtx itsEnergy;
  Mtx itsStartValues;

  MultivarFunction& itsObjective;

public:
  AnnealingOptimizer(MultivarFunction& objective, AnnealOpts& opts)
    :
    itsOpts(opts),
    itsRunNum(0),
    itsDeltas(itsOpts.deltas),
    itsNumFunEvals(itsOpts.numRuns, 1),
    itsBestCosts(itsOpts.numRuns, 1),
    itsModelHist(itsOpts.numModelParams, int(itsOpts.coolingSchedule.sum())),
    itsBestModels(itsOpts.numModelParams, itsOpts.numRuns),
    itsEnergy(int(itsOpts.coolingSchedule.sum()), itsOpts.numRuns),
    itsStartValues(itsOpts.centers),

    itsObjective(objective)
  {
    itsEnergy.setAll(std::numeric_limits<double>::max());
  }

  void doOneRun();

  void optimize()
  {
    for (int i = 0; i < itsOpts.numRuns; ++i)
      doOneRun();
  }

  mxArray* getOutput()
  {
    const char* fieldnames[] =
      {
        "bestCost",
        "mhat",
        "model",
        "energy",
        "numFunEvals"
      };

    mxArray* output = mxCreateStructMatrix(1, 1, 5, fieldnames);

    int best_pos = -1;
    double best_cost = itsBestCosts.min(&best_pos);

    Mtx mhat = itsBestModels.column(best_pos);

    mxSetField(output, 0, "bestCost", mxCreateScalarDouble(best_cost));
    mxSetField(output, 0, "mhat", mhat.makeMxArray());
    mxSetField(output, 0, "model", itsModelHist.makeMxArray());
    mxSetField(output, 0, "energy", itsEnergy.makeMxArray());
    mxSetField(output, 0, "numFunEvals", itsNumFunEvals.makeMxArray());

    return output;
  }

  // Make a new set of models, based on srcModel, but with a range of new
  // parameters values for whichParam
  Mtx makeTestModels(int whichParam, const Mtx& srcModel)
  {
    DOTRACE("makeTestModels");

    const double delta = itsDeltas.at(whichParam);

    const double current_param_val = srcModel.at(whichParam);

    const double lowerbound = itsOpts.bounds.at(whichParam,0);
    const double upperbound = itsOpts.bounds.at(whichParam,1);

    int num_testmodels = 0;
    {
      for (MtxConstIter iter = itsOpts.valueScalingRange.rowIter(0);
           iter.hasMore(); ++iter)
        {
          double value = current_param_val + (*iter) * delta;
          if ( (value > lowerbound) && (value < upperbound) )
            ++num_testmodels;
        }
    }

    Mtx newModels(srcModel.mrows(), num_testmodels);

    const Slice srcModelCol = srcModel.column(0);

    for (int i = 0; i < num_testmodels; ++i)
      {
        newModels.column(i) = srcModelCol;
      }

    {
      int col = 0;
      for (MtxConstIter iter = itsOpts.valueScalingRange.rowIter(0);
           iter.hasMore(); ++iter)
        {
          double value = current_param_val + (*iter) * delta;
          if ( (value > lowerbound) && (value < upperbound) )
            newModels.at(whichParam, col++) = value;
        }
    }

    return newModels;
  }


  // Returns the cost at the best model that was found
  double visitParameters(Mtx& bestModel, const double temp);

  Mtx makeRandomModels()
  {
    Mtx models = matlabRand(itsOpts.numModelParams,
                            itsOpts.numStartingPoints);

    models -= 0.5;
    models *= 2;

    for (int r = 0; r < models.mrows(); ++r)
      {
        models.row(r) *= itsDeltas.at(r);
        models.row(r) += itsStartValues.at(r);
      }

    return models;
  }

  Mtx createStartingModel(double* criticalTemp)
  {
    DOTRACE("createStartingModel");

    const Mtx startingModels = makeRandomModels();

    const Mtx startingCosts = itsObjective.evaluateEach(startingModels);

    *criticalTemp =
      log10(startingCosts.sum() / startingCosts.nelems())
      - itsOpts.tempScales.at(itsRunNum);

    int startingPos = 0;
    const double startingCost = startingCosts.min(&startingPos);

    if (itsOpts.talking)
      {
        std::cerr << "\nStarting cost "
                  << std::setw(7) << std::fixed << std::setprecision(2)
                  << startingCost;
        std::cerr << "\n\nBeginning run #" << itsRunNum+1 << ". ";
        std::cerr << "Critical temperature at "
                  << std::setw(3) << std::fixed << std::setprecision(2)
                  << pow(10.0, *criticalTemp) << ".\n";

        std::cerr << "------------------------------------------------\n\n";
        std::cerr << "f-Calls\t\tTemperature\tMinimum f-Value\n";
        std::cerr << "------------------------------------------------\n";
      }

    return Mtx(startingModels.column(startingPos));
  }

  void displayParams(const Mtx& model, double cost)
  {
    if (!itsOpts.talking) return;

    std::cerr << "\nparams: ";
    for (int i = 0; i < model.nelems(); ++i)
      std::cerr << std::fixed << std::setprecision(6) << model.at(i) << " ";
    std::cerr << "\ncost: "
              << std::setw(7) << std::fixed << std::setprecision(4)
              << cost
              << "\n";
  }

  void runSimplex()
  {
    SimplexOptimizer opt(itsObjective,
                         Mtx(itsBestModels.column(itsRunNum)),
                         fstring("notify"),
                         itsOpts.numModelParams,
                         10000000, // maxFunEvals
                         10000, // maxIter
                         1e-4, // tolx
                         1e-4 // tolf
                         );

    /*int exitFlag =*/ opt.optimize();

    double Ostar = opt.bestFval();

    Mtx mstar = opt.bestParams();

    if (Ostar < itsBestCosts.at(itsRunNum))
      {
        displayParams(mstar, Ostar);

        if (itsOpts.talking)
          std::cerr << opt.iterCount() << " iterations\n";

        bool inBounds = true;

        for (int i = 0; i < mstar.nelems(); ++i)
          {
            double val = mstar.at(i);
            if (val < itsOpts.bounds.at(i, 0)) { inBounds = false; break; }
            if (val > itsOpts.bounds.at(i, 1)) { inBounds = false; break; }
          }

        if (inBounds)
          {
            itsBestModels.column(itsRunNum) = mstar;
            itsBestCosts.at(itsRunNum) = Ostar;
            if (itsOpts.talking)
              std::cerr << "\nSimplex method lowered cost "
                        << "and remained within constraints.\n\n";
          }
        else
          {
            if (itsOpts.talking)
              std::cerr << "\nSimplex method lowered cost "
                        << "but failed to remain within constraints.\n\n";
          }
      }
  }
};

double AnnealingOptimizer::visitParameters(Mtx& bestModel, const double temp)
{
  int nevals = 0;
  int s = -1;

  Mtx costs(0,0);

  for (int x = 0; x < itsDeltas.nelems(); ++x)
    {
      const double delta = itsDeltas.at(x);

      if (delta == 0.0) continue;

      Mtx modelmatrix = makeTestModels(x, bestModel);

      costs = itsObjective.evaluateEach(modelmatrix);

      nevals += costs.nelems();

      // Sample from probability distribution
      s = sampleFromPdf(temp, costs);

      bestModel.at(x, 0) = modelmatrix.at(x, s);
    }

  if (s < 0)
    {
      throw Util::Error("didn't run any annealing iterations because "
                        "there were no non-zero deltas");
    }

  itsNumFunEvals.at(itsRunNum) += nevals;

  return costs.at(s);
}

//---------------------------------------------------------------------
//
// AnnealingOptimizer::doOneRun()
//
//---------------------------------------------------------------------

void AnnealingOptimizer::doOneRun()
{
DOTRACE("AnnealingOptimizer::doOneRun");

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

  int nvisits = 0;
  double criticalTemp = 0.0;

  Mtx bestModel = createStartingModel(&criticalTemp);

  Mtx minUsedParams = bestModel;
  Mtx maxUsedParams = bestModel;

  for (int temps_i = 0; temps_i < itsOpts.numTemps; ++temps_i)
    {
      // so the temperatures range from 10^(crit_temp+1) ... 10^(crit_temp-1)
      const double temp =
        pow(10.0, (criticalTemp + 1.0
                   - 2.0*double(temps_i)/(itsOpts.numTemps-1)));

      const int temp_repeat = int(itsOpts.coolingSchedule.at(temps_i));

      for (int repeat = 0; repeat < temp_repeat; ++repeat)
        {
          ++nvisits;

          if (itsOpts.talking && (nvisits % 10 == 0))
            {
              std::cerr << std::setw(7) << int(itsNumFunEvals.at(itsRunNum))
                        << "\t\t"
                        << std::setw(7) << std::fixed << std::setprecision(2)
                        << temp
                        << "\t\t"
                        << std::setw(7) << std::fixed << std::setprecision(2)
                        << itsEnergy.column(itsRunNum).min()
                        << "\n";
            }

          itsEnergy.at(nvisits-1,itsRunNum) =
            visitParameters(bestModel, temp);

          for (int i = 0; i < bestModel.nelems(); ++i)
            {
              minUsedParams.at(i) = std::min(double(minUsedParams.at(i)),
                                             double(bestModel.at(i)));

              maxUsedParams.at(i) = std::max(double(maxUsedParams.at(i)),
                                             double(bestModel.at(i)));
            }

          itsModelHist.column(nvisits-1) = bestModel;
        }
    }


  // Update deltas
  for (int i = 0; i < itsDeltas.nelems(); ++i)
    {
      itsDeltas.at(i) = 0.75 * (maxUsedParams.at(i) - minUsedParams.at(i));
    }

  // Update bests FIXME ought to use a smarter algorithm to keep track of the
  // best cost

  Mtx currentEnergy = itsEnergy.column(itsRunNum);

  int best_pos = 0;
  const double best_energy = currentEnergy.min(&best_pos);

  itsBestCosts.at(itsRunNum) = best_energy;
  itsBestModels.column(itsRunNum) = itsModelHist.column(best_pos);

  displayParams(itsBestModels.column(itsRunNum), itsBestCosts.at(itsRunNum));

  if (itsOpts.doNewton) runSimplex();

  itsStartValues = itsBestModels.column(itsRunNum);

  ++itsRunNum;
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
      throw Util::Error("Error: annealHelper was called with more "
                        "than the declared number of outputs (1).");
    }

  const int NDECLARED = 2;

  if (nrhs < NDECLARED)
    {
      throw Util::Error("Error: annealHelper was called with fewer "
                        "than the declared number of inputs (7).");
    }

  mlfEnterNewContext(0, NDECLARED, prhs[0], prhs[1]);

  int nvararg = nrhs - NDECLARED;
  mxArray** pvararg = prhs + NDECLARED;

  try
    {
      AnnealOpts opts(prhs[0]);

      MatlabFunction objective(MxWrapper::extractString(prhs[1]), // funcName
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
