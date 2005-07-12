///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 2001-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Thu Mar  8 09:34:12 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

#include "cmodel/classifier.h"

#include "mtx/mathspecial.h"

#include "optim/annealingoptimizer.h"
#include "optim/multivarfunction.h"
#include "optim/simplexoptimizer.h"

#include "rutz/error.h"
#include "rutz/fstring.h"

#include <limits>
#include <matrix.h>

#include "rutz/trace.h"
#include "rutz/debug.h"

using rutz::fstring;

class LLEvaluator : public MultivarFunction
{
  Classifier& itsC;
  const mtx& itsObjs;
  const mtx& itsInc;

protected:
  virtual double doEvaluate(const mtx& x)
  {
    GVX_TRACE("LLEvaluator::doEvaluate");
    slice params(x.column(0));
    return -itsC.currentLogL(params, itsObjs, itsInc);
  }

public:
  LLEvaluator(Classifier& c, const mtx& objs, const mtx& inc) :
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
  const mtx& testSize(const mtx& objParams)
  {
    GVX_TRACE("<classifier.cc>::testSize");
    if (objParams.ncols() != Classifier::DIM_OBJ_PARAMS+1)
      throw rutz::error(fstring("objParams must have "
                                "DIM_OBJ_PARAMS+1 columns "
                                "(expected ", Classifier::DIM_OBJ_PARAMS+1,
                                ", got ", objParams.ncols(), ")"), SRC_POS);

    return objParams;
  }
}

Classifier::Classifier(const mtx& objParams) :
  itsObjParams(testSize(objParams)),
  itsObjCategories(objParams(col_range_n(0, 1))),
  itsObjects(objParams(col_range_n(1, DIM_OBJ_PARAMS))),
  itsNumAllExemplars(objParams.mrows()),
  itsDiffEvidence(mtx::zeros(itsNumAllExemplars,1)),
  itsObservedIncidenceCache(mtx::empty_mtx()),
  itsCachedLogL_1_2(0.0)
{
GVX_TRACE("Classifier::Classifier");
}

Classifier::~Classifier()
{
GVX_TRACE("Classifier::~Classifier");
}

int Classifier::numModelParams() const
{
GVX_TRACE("Classifier::numModelParams");

  return DIM_OBJ_PARAMS // one attentional weight per obj param
    + 1 // plus one param for the threshold
    + 1; // plus one param for the std deviation of the noise
}

mtx Classifier::modelParamsBounds() const
{
GVX_TRACE("Classifier::modelParamsBounds");

  mtx bounds = mtx::zeros(numModelParams(), 2);

  const int row = fillModelParamsBounds(bounds, 0);

  if (row != numModelParams())
    {
      throw rutz::error(fstring("not all rows were filled "
                                "in modelParamsBounds() "
                                "(expected ", numModelParams(),
                                ", got ", row, ")"), SRC_POS);
    }

  return bounds;
}

mtx Classifier::forwardProbit(const mtx& diffEv,
                              double thresh, double sigmaNoise)
{
GVX_TRACE("Classifier::forwardProbit");

  const double divisor = (1.0 / M_SQRT2) * (1.0 / sigmaNoise);

  mtx_const_iter diffev = diffEv.column_iter(0);

  mtx pp = mtx::zeros(diffEv.mrows(), 1);
  mtx_iter ppiter = pp.column_iter(0);

  // alpha = (thresh - diffEvidence) / sigmaNoise
  //
  // pp = 0.5 * erfc(alpha / sqrt(2))

  for (; ppiter.has_more(); ++diffev, ++ppiter) {
    double alpha_val = thresh - *diffev;

    *ppiter = 0.5*dash::erfc(alpha_val * divisor);
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

double Classifier::computeLogL(const mtx& predictedProbability,
                               const mtx& observedIncidence)
{
GVX_TRACE("Classifier::computeLogL");

  if (observedIncidence != itsObservedIncidenceCache ||
      itsCachedLogL_1_2 == 0.0)
    {
      itsCachedLogL_1_2 = 0.0;

      for(int k = 0; k < observedIncidence.mrows(); ++k)
        {
          double oi1 = observedIncidence.at(k,0);
          double oi2 = observedIncidence.at(k,1);

          // term 1
          itsCachedLogL_1_2 += dash::gammaln(1.0 + oi1 + oi2);

          // term 2
          itsCachedLogL_1_2 -= dash::gammaln(1.0+oi1);
          itsCachedLogL_1_2 -= dash::gammaln(1.0+oi2);
        }

      itsObservedIncidenceCache = observedIncidence;
    }

  double logL_3 = 0.0;

  const double LOG_10_MINUS_50 = -115.1293;

  mtx_const_iter oi1iter = observedIncidence.column_iter(0);
  mtx_const_iter oi2iter = observedIncidence.column_iter(1);

  mtx_const_iter ppiter = predictedProbability.column_iter(0);

  for(; ppiter.has_more(); ++ppiter, ++oi1iter, ++oi2iter)
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
mtx Classifier::classifyObjects(slice& modelParams, const mtx& testObjects)
{
GVX_TRACE("Classifier::classifyObjects");

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

double Classifier::currentLogL(slice& modelParams,
                               const mtx& testObjects,
                               const mtx& observedIncidence)
{
GVX_TRACE("Classifier::currentLogL");

  mtx pp = classifyObjects(modelParams, testObjects);

  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.

  return computeLogL(pp, observedIncidence);
}

double Classifier::fullLogL(const mtx& observedIncidence)
{
GVX_TRACE("Classifier::fullLogL");

  mtx observedProb = mtx::zeros(numAllExemplars(), 1);
  mtx_iter opiter = observedProb.column_iter(0);

  mtx_const_iter oi1iter = observedIncidence.column_iter(0);
  mtx_const_iter oi2iter = observedIncidence.column_iter(1);

  for (; opiter.has_more(); ++opiter, ++oi1iter, ++oi2iter)
    *opiter = (*oi1iter / (*oi1iter + *oi2iter));

  return computeLogL(observedProb, observedIncidence);
}

double Classifier::deviance(slice& modelParams,
                            const mtx& testObjects,
                            const mtx& observedIncidence)
{
GVX_TRACE("Classifier::deviance");

  double llc = currentLogL(modelParams, testObjects, observedIncidence);
  double llf = fullLogL(observedIncidence);

  return -2 * (llc - llf);
}

namespace
{
  mtx getTestObjects(const mx_wrapper& extraArgs)
  {
    mtx result = extraArgs.get_mtx_field("testObjects");

    if (result.ncols() != Classifier::DIM_OBJ_PARAMS)
      {
        throw rutz::error("wrong number of columns in 'testObjects' "
                          "(should be equal to DIM_OBJ_PARAMS)", SRC_POS);
      }

    return result;
  }

  void checkModelParams(const mtx& allModelParams,
                        const int expectedNumber)
  {
    if (allModelParams.mrows() != expectedNumber)
      {
        throw rutz::error(fstring("wrong number of model params "
                                  "(expected ", expectedNumber,
                                  ", got ", allModelParams.mrows(), ")"),
                          SRC_POS);
      }
  }
}

Classifier::RequestResult
Classifier::handleRequest(rutz::fstring action,
                          const mtx& allModelParams,
                          const mx_wrapper& extraArgs)
{
GVX_TRACE("Classifier::handleRequest");

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
      GVX_TRACE("Classifier::handleRequest-bounds");

      return modelParamsBounds();
    }

  else if ( action == "ll" || action == "llc" )
    {
      GVX_TRACE("Classifier::handleRequest-llc");

      const mtx observedIncidence = extraArgs.get_mtx_field("observedIncidence");

      const mtx testObjects = getTestObjects(extraArgs);

      checkModelParams(allModelParams, numModelParams());

      mtx result = mtx::zeros(allModelParams.ncols(), 1);
      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          slice modelParams(allModelParams.column(i));
          result.at(i) = multiplier *
            currentLogL(modelParams, testObjects, observedIncidence);
        }

      return result;
    }

  else if ( action == "llf" )
    {
      GVX_TRACE("Classifier::handleRequest-llf");

      const mtx observedIncidence = extraArgs.get_mtx_field("observedIncidence");

      checkModelParams(allModelParams, numModelParams());

      mtx result = mtx::zeros(allModelParams.ncols(), 1);
      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          result.at(i) = multiplier * fullLogL(observedIncidence);
        }

      return result;
    }

  else if ( action == "dev" )
    {
      GVX_TRACE("Classifier::handleRequest-dev");

      const mtx observedIncidence = extraArgs.get_mtx_field("observedIncidence");

      const mtx testObjects = getTestObjects(extraArgs);

      checkModelParams(allModelParams, numModelParams());

      mtx result = mtx::zeros(allModelParams.ncols(), 1);
      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          slice modelParams(allModelParams.column(i));
          result.at(i) = multiplier *
            deviance(modelParams, testObjects, observedIncidence);
        }

      return result;
    }

  else if ( action == "classify" )
    {
      GVX_TRACE("Classifier::handleRequest-classify");

      const mtx testObjects = getTestObjects(extraArgs);

      checkModelParams(allModelParams, numModelParams());

      mtx result = mtx::zeros(testObjects.mrows(), allModelParams.ncols());

      for (int i = 0; i < allModelParams.ncols(); ++i)
        {
          slice modelParams(allModelParams.column(i));
          result.column(i) = this->classifyObjects(modelParams, testObjects);
        }

      return result;
    }

  else if ( action == "simplex" )
    {
      GVX_TRACE("Classifier::handleRequest-simplex");

      const mtx observedIncidence = extraArgs.get_mtx_field("observedIncidence");

      const mtx testObjects = getTestObjects(extraArgs);

      LLEvaluator objective(*this, testObjects, observedIncidence);

      checkModelParams(allModelParams, numModelParams());

      SimplexOptimizer opt(objective,
                           allModelParams.as_column(),
                           "notify",
                           allModelParams.nelems(),
                           extraArgs.get_int_field("maxfun"),
                           extraArgs.get_int_field("maxiter")
                           );

      int exitFlag = opt.optimize();

      mx_wrapper result;

      result.set_field("bestParams", opt.bestParams());
      result.set_field("bestFval", opt.bestFval());
      result.set_field("exitFlag", exitFlag);
      result.set_field("iterations", opt.iterCount());
      result.set_field("funcCount", opt.funcCount());
      result.set_field("algorithm", opt.algorithm());

      return result;
    }

  else if ( action == "anneal" )
    {
      GVX_TRACE("Classifier::handleRequest-anneal");

      const mtx observedIncidence = extraArgs.get_mtx_field("observedIncidence");

      const mtx testObjects = getTestObjects(extraArgs);

      LLEvaluator objective(*this, testObjects, observedIncidence);

      AnnealOpts opts(extraArgs.get_field("annealOpts"));

      AnnealingOptimizer ar(objective, opts);

      ar.optimize();

      mx_wrapper result;
      ar.getOutput(result);
      return result;
    }

  return RequestResult();
}

// Count the category training exemplars
int Classifier::countCategory(int category) const
{
GVX_TRACE("Classifier::countCategory");
  int n = 0;
  mtx_const_iter iter = itsObjCategories.column_iter(0);
  for (; iter.has_more(); ++iter)
    {
      if (*iter == category)
        ++n;
    }
  return n;
}

mtx Classifier::objectsOfCategory(int category) const
{
GVX_TRACE("Classifier::objectsOfCategory");

  const int nobjs = countCategory(category);

  if (nobjs == 0)
    throw rutz::error(fstring("no objects found in category ", category),
                      SRC_POS);

  mtx result = mtx::zeros(nobjs, DIM_OBJ_PARAMS);

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
GVX_TRACE("Classifier::exemplarCategory");
  return int(itsObjCategories.at(i));
}

slice Classifier::exemplar(int i) const
{
GVX_TRACE("Classifier::exemplar");
  return itsObjects.row(i);
}

int Classifier::fillModelParamsBounds(mtx& bounds, int startRow) const
{
GVX_TRACE("Classifier::fillModelParamsBounds");

  const double minus_inf = -std::numeric_limits<double>::max();
  const double plus_inf = std::numeric_limits<double>::max();

  // All lower bounds default to "minus infinity"
  // All upper bounds default to "plus infinity"
  bounds.sub(col_range_n(0,1)).clear(minus_inf);
  bounds.sub(col_range_n(1,1)).clear(plus_inf);

  return DIM_OBJ_PARAMS+2;
}

static const char vcid_classifier_cc[] = "$Id$ $URL$";
#endif // !CLASSIFIER_CC_DEFINED
