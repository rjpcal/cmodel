///////////////////////////////////////////////////////////////////////
//
// annealingoptimizer.h
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Tue Feb 19 09:59:33 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef ANNEALINGOPTIMIZER_H_DEFINED
#define ANNEALINGOPTIMIZER_H_DEFINED

#include "mtx/mtx.h"

class MultivarFunction;

///////////////////////////////////////////////////////////////////////
//
// AnnealOpts
//
///////////////////////////////////////////////////////////////////////

struct AnnealOpts
{
private:
  Mtx makeScalingRange(int gridSpacing);

  Mtx makeCoolingSchedule(int scale);

public:
  AnnealOpts(const mxArray* arr);

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
  AnnealingOptimizer(MultivarFunction& objective, AnnealOpts& opts);

  void doOneRun();

  void optimize()
  {
    for (int i = 0; i < itsOpts.numRuns; ++i)
      doOneRun();
  }

  mxArray* getOutput();

  // Make a new set of models, based on srcModel, but with a range of new
  // parameters values for whichParam
  Mtx makeTestModels(int whichParam, const Mtx& srcModel);

  // Returns the cost at the best model that was found
  double visitParameters(Mtx& bestModel, const double temp);

  Mtx makeRandomModels();

  Mtx createStartingModel(double* criticalTemp);

  void displayParams(const Mtx& model, double cost);

  void runSimplex();
};

static const char vcid_annealingoptimizer_h[] = "$Header$";
#endif // !ANNEALINGOPTIMIZER_H_DEFINED
