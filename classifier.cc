///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Thu Mar  8 09:34:12 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

#include "cmodel/classifier.h"

#include "mtx/num.h"

#include "optim/annealingoptimizer.h"
#include "optim/multivarfunction.h"
#include "optim/simplexoptimizer.h"

#include "util/error.h"
#include "util/strings.h"

#include <limits>
#include <matrix.h>

#include "util/trace.h"
#include "util/debug.h"

class LLEvaluator : public MultivarFunction
{
  Classifier& itsC;
  const Mtx& itsObjs;
  const Mtx& itsInc;

protected:
  virtual double doEvaluate(const Mtx& x)
  {
    DOTRACE("LLEvaluator::doEvaluate");
    Slice params(x.column(0));
    return -itsC.currentLogL(params, itsObjs, itsInc);
  }

public:
  LLEvaluator(Classifier& c, const Mtx& objs, const Mtx& inc) :
    itsC(c), itsObjs(objs), itsInc(inc) {}
};

///////////////////////////////////////////////////////////////////////
//
// Classifier member definitions
//
///////////////////////////////////////////////////////////////////////

namespace
{
  // Just check that objParams has the right size
  const Mtx& testSize(const Mtx& objParams)
  {
    DOTRACE("<classifier.cc>::testSize");
    if (objParams.ncols() != Classifier::DIM_OBJ_PARAMS+1)
      throw Util::Error(fstring("objParams must have "
                                "DIM_OBJ_PARAMS+1 columns "
                                "(expected ", Classifier::DIM_OBJ_PARAMS+1,
                                ", got ", objParams.ncols(), ")"));

    return objParams;
  }
}

Classifier::Classifier(const Mtx& objParams) :
  itsObjParams(testSize(objParams)),
  itsObjCategories(objParams(col_range_n(0, 1))),
  itsObjects(objParams(col_range_n(1, DIM_OBJ_PARAMS))),
  itsNumAllExemplars(objParams.mrows()),
  itsDiffEvidence(itsNumAllExemplars,1),
  itsObservedIncidenceCache(0,0),
  itsCachedLogL_1_2(0.0)
{
DOTRACE("Classifier::Classifier");
}

Classifier::~Classifier()
{
DOTRACE("Classifier::~Classifier");
}

int Classifier::numModelParams() const
{
DOTRACE("Classifier::numModelParams");

  return DIM_OBJ_PARAMS // one attentional weight per obj param
    + 1 // plus one param for the threshold
    + 1; // plus one param for the std deviation of the noise
}

Mtx Classifier::modelParamsBounds() const
{
DOTRACE("Classifier::modelParamsBounds");

  Mtx bounds(numModelParams(), 2);

  const int row = fillModelParamsBounds(bounds, 0);

  if (row != numModelParams())
    {
      throw Util::Error(fstring("not all rows were filled "
                                "in modelParamsBounds() "
                                "(expected ", numModelParams(),
                                ", got ", row, ")"));
    }

  return bounds;
}

Mtx Classifier::forwardProbit(const Mtx& diffEv,
                              double thresh, double sigmaNoise)
{
DOTRACE("Classifier::forwardProbit");

  const double divisor = (1.0 / Num::SQRT_2) * (1.0 / sigmaNoise);

  MtxConstIter diffev = diffEv.columnIter(0);

  Mtx pp(diffEv.mrows(), 1);
  MtxIter ppiter = pp.columnIter(0);

  // alpha = (thresh - diffEvidence) / sigmaNoise
  //
  // pp = 0.5 * erfc(alpha / sqrt(2))

  for (; ppiter.hasMore(); ++diffev, ++ppiter) {
    double alpha_val = thresh - *diffev;

    *ppiter = 0.5*Num::erfc(alpha_val * divisor);
  }

  return pp;
}

//---------------------------------------------------------------------
//
// ll = sum(gammaln(1+sum(observedIncidence,2))) - ...
//      sum(sum(gammaln(observedIncidence+1)))  + ...
//      sum(sum(observedIncidence.*log(predictedProbability)));
//
//---------------------------------------------------------------------

double Classifier::computeLogL(const Mtx& predictedProbability,
                               const Mtx& observedIncidence)
{
DOTRACE("Classifier::computeLogL");

  if (observedIncidence != itsObservedIncidenceCache ||
      itsCachedLogL_1_2 == 0.0)
    {
      itsCachedLogL_1_2 = 0.0;

      for(int k = 0; k < observedIncidence.mrows(); ++k)
        {
          double oi1 = observedIncidence.at(k,0);
          double oi2 = observedIncidence.at(k,1);

          // term 1
          itsCachedLogL_1_2 += Num::gammaln(1.0 + oi1 + oi2);

          // term 2
          itsCachedLogL_1_2 -= Num::gammaln(1.0+oi1);
          itsCachedLogL_1_2 -= Num::gammaln(1.0+oi2);
        }

      itsObservedIncidenceCache = observedIncidence;
    }

  double logL_3 = 0.0;

  const double LOG_10_MINUS_50 = -115.1293;

  MtxConstIter oi1iter = observedIncidence.columnIter(0);
  MtxConstIter oi2iter = observedIncidence.columnIter(1);

  MtxConstIter ppiter = predictedProbability.columnIter(0);

  for(; ppiter.hasMore(); ++ppiter, ++oi1iter, ++oi2iter)
    {
      double oi1 = *oi1iter;
      double oi2 = *oi2iter;

      // term3
      double pp_val = *ppiter;

      if (pp_val < 1e-50) logL_3 += oi1 * LOG_10_MINUS_50;
      else                logL_3 += oi1 * log(pp_val);

      pp_val = 1.0 - pp_val;

      if (pp_val < 1e-50) logL_3 += oi2 * LOG_10_MINUS_50;
      else                logL_3 += oi2 * log(pp_val);
    }

  return itsCachedLogL_1_2 + logL_3;
}

// Returns the classification probability for each of 'objects'
Mtx Classifier::classifyObjects(Slice& modelParams, const Mtx& testObjects)
{
DOTRACE("Classifier::classifyObjects");

  itsDiffEvidence.resize(testObjects.mrows(), 1);

  computeDiffEv(testObjects, modelParams, itsDiffEvidence);

  //---------------------------------------------------------------------
  //
  // Compute the predicted probabilities based on the similarity
  // differences in diffEvidence.
  //

  // the threshold is right after the DIM_OBJ_PARAMS number of
  // attentional weights
  const double thresh = modelParams[DIM_OBJ_PARAMS];

  const double sigmaNoise = (modelParams.nelems() > DIM_OBJ_PARAMS+1) ?
    computeSigmaNoise(modelParams[DIM_OBJ_PARAMS+1]) : 1.0;

  // predictedProbability = forwardProbit(diffEvidence, thresh, sigmaNoise);
  return forwardProbit(itsDiffEvidence, thresh, sigmaNoise);
}

double Classifier::currentLogL(Slice& modelParams,
                               const Mtx& testObjects,
                               const Mtx& observedIncidence)
{
DOTRACE("Classifier::currentLogL");

  Mtx pp = classifyObjects(modelParams, testObjects);

  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.

  return computeLogL(pp, observedIncidence);
}

double Classifier::fullLogL(const Mtx& observedIncidence)
{
DOTRACE("Classifier::fullLogL");

  Mtx observedProb(numAllExemplars(), 1);
  MtxIter opiter = observedProb.columnIter(0);

  MtxConstIter oi1iter = observedIncidence.columnIter(0);
  MtxConstIter oi2iter = observedIncidence.columnIter(1);

  for (; opiter.hasMore(); ++opiter, ++oi1iter, ++oi2iter)
    *opiter = (*oi1iter / (*oi1iter + *oi2iter));

  return computeLogL(observedProb, observedIncidence);
}

double Classifier::deviance(Slice& modelParams,
                            const Mtx& testObjects,
                            const Mtx& observedIncidence)
{
DOTRACE("Classifier::deviance");

  double llc = currentLogL(modelParams, testObjects, observedIncidence);
  double llf = fullLogL(observedIncidence);

  return -2 * (llc - llf);
}

namespace
{
  Mtx getTestObjects(const MxWrapper& extraArgs)
  {
    Mtx result = extraArgs.getField("testObjects").asMtx();

    if (result.ncols() != Classifier::DIM_OBJ_PARAMS)
      {
        throw Util::Error("wrong number of columns in 'testObjects' "
                          "(should be equal to DIM_OBJ_PARAMS)");
      }

    return result;
  }

  void checkModelParams(const Mtx& allModelParams,
                        const int expectedNumber)
  {
    if (allModelParams.mrows() != expectedNumber)
      {
        throw Util::Error(fstring("wrong number of model params "
                                  "(expected ", expectedNumber,
                                  ", got ", allModelParams.mrows(), ")"));
      }
  }
}

Classifier::RequestResult
Classifier::handleRequest(fstring action,
                          const Mtx& allModelParams,
                          const MxWrapper& extraArgs)
{
DOTRACE("Classifier::handleRequest");

  int multiplier = 1;
  // check for minus sign
  if (action.c_str()[0] == '-')
    {
      multiplier = -1;
      fstring trimmed = action.c_str() + 1;
      action = trimmed;
    }

  //---------------------------------------------------------------------
  //
  // Call the requested computational function for each set of model
  // params
  //
  //---------------------------------------------------------------------

  if ( action == "bounds" )
    {
      DOTRACE("Classifier::handleRequest-bounds");

      return modelParamsBounds();
    }

  else if ( action == "ll" || action == "llc" )
    {
      DOTRACE("Classifier::handleRequest-llc");

      const Mtx observedIncidence = extraArgs.getField("observedIncidence").asMtx();

      const Mtx testObjects = getTestObjects(extraArgs);

      checkModelParams(allModelParams, numModelParams());

      Mtx result(allModelParams.ncols(), 1);
      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          Slice modelParams(allModelParams.column(i));
          result.at(i) = multiplier *
            currentLogL(modelParams, testObjects, observedIncidence);
        }

      return result;
    }

  else if ( action == "llf" )
    {
      DOTRACE("Classifier::handleRequest-llf");

      const Mtx observedIncidence = extraArgs.getField("observedIncidence").asMtx();

      checkModelParams(allModelParams, numModelParams());

      Mtx result(allModelParams.ncols(), 1);
      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          result.at(i) = multiplier * fullLogL(observedIncidence);
        }

      return result;
    }

  else if ( action == "dev" )
    {
      DOTRACE("Classifier::handleRequest-dev");

      const Mtx observedIncidence = extraArgs.getField("observedIncidence").asMtx();

      const Mtx testObjects = getTestObjects(extraArgs);

      checkModelParams(allModelParams, numModelParams());

      Mtx result(allModelParams.ncols(), 1);
      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          Slice modelParams(allModelParams.column(i));
          result.at(i) = multiplier *
            deviance(modelParams, testObjects, observedIncidence);
        }

      return result;
    }

  else if ( action == "classify" )
    {
      DOTRACE("Classifier::handleRequest-classify");

      const Mtx testObjects = getTestObjects(extraArgs);

      checkModelParams(allModelParams, numModelParams());

      Mtx result(testObjects.mrows(), allModelParams.ncols());

      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          Slice modelParams(allModelParams.column(i));
          result.column(i) = this->classifyObjects(modelParams, testObjects);
        }

      return result;
    }

  else if ( action == "simplex" )
    {
      DOTRACE("Classifier::handleRequest-simplex");

      const Mtx observedIncidence = extraArgs.getField("observedIncidence").asMtx();

      const Mtx testObjects = getTestObjects(extraArgs);

      LLEvaluator objective(*this, testObjects, observedIncidence);

      checkModelParams(allModelParams, numModelParams());

      SimplexOptimizer opt(objective,
                           allModelParams.asColumn(),
                           "notify",
                           allModelParams.nelems(),
                           extraArgs.getField("maxfun").asInt(),
                           extraArgs.getField("maxiter").asInt()
                           );

      int exitFlag = opt.optimize();

      MxWrapper result;

      result.setField("bestParams", opt.bestParams());
      result.setField("bestFval", opt.bestFval());
      result.setField("exitFlag", exitFlag);
      result.setField("iterations", opt.iterCount());
      result.setField("funcCount", opt.funcCount());
      result.setField("algorithm", opt.algorithm());

      return result;
    }

  else if ( action == "anneal" )
    {
      DOTRACE("Classifier::handleRequest-anneal");

      const Mtx observedIncidence = extraArgs.getField("observedIncidence").asMtx();

      const Mtx testObjects = getTestObjects(extraArgs);

      LLEvaluator objective(*this, testObjects, observedIncidence);

      MxWrapper annealArgs(extraArgs.getField("annealOpts"));

      // FIXME this is ugly
      mxArray* args = annealArgs.release();
      AnnealOpts opts(args);
      mxDestroyArray(args);

      AnnealingOptimizer ar(objective, opts);

      ar.optimize();

      // FIXME this is also ugly
      mxArray* output = ar.getOutput();
      MxWrapper result(output);
      mxDestroyArray(output);
      return result;
    }

  return RequestResult();
}

// Count the category training exemplars
int Classifier::countCategory(int category) const
{
DOTRACE("Classifier::countCategory");
  int n = 0;
  MtxConstIter iter = itsObjCategories.columnIter(0);
  for (; iter.hasMore(); ++iter)
    {
      if (*iter == category)
        ++n;
    }
  return n;
}

Mtx Classifier::objectsOfCategory(int category) const
{
DOTRACE("Classifier::objectsOfCategory");

  const int nobjs = countCategory(category);

  if (nobjs == 0)
    throw Util::Error(fstring("no objects found in category ", category));

  Mtx result(nobjs, DIM_OBJ_PARAMS);

  int r = 0;
  for (int i = 0; i < itsObjects.mrows(); ++i)
    {
      if (exemplarCategory(i) == category)
        result.row(r++) = exemplar(i);
    }

  return result;
}

int Classifier::exemplarCategory(int i) const
{
DOTRACE("Classifier::exemplarCategory");
  return int(itsObjCategories.at(i));
}

Slice Classifier::exemplar(int i) const
{
DOTRACE("Classifier::exemplar");
  return itsObjects.row(i);
}

int Classifier::fillModelParamsBounds(Mtx& bounds, int startRow) const
{
DOTRACE("Classifier::fillModelParamsBounds");

  const double minus_inf = -std::numeric_limits<double>::max();
  const double plus_inf = std::numeric_limits<double>::max();

  // All lower bounds default to "minus infinity"
  // All upper bounds default to "plus infinity"
  bounds.sub(col_range_n(0,1)).setAll(minus_inf);
  bounds.sub(col_range_n(1,1)).setAll(plus_inf);

  return DIM_OBJ_PARAMS+2;
}

static const char vcid_classifier_cc[] = "$Header$";
#endif // !CLASSIFIER_CC_DEFINED
