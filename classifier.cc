///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:34:12 2001
// written: Fri Apr  6 10:51:42 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

#include "classifier.h"

#include "error.h"
#include "num.h"
#include "strings.h"

#include "trace.h"

///////////////////////////////////////////////////////////////////////
//
// Classifier member definitions
//
///////////////////////////////////////////////////////////////////////

Classifier::Classifier(const Mtx& objParams,
							  const Mtx& observedIncidence) :
  itsObjCategories(objParams.columns(0, 1)),
  itsObjects(objParams.columns(1, DIM_OBJ_PARAMS)),
  itsNumAllExemplars(objParams.mrows()),
  itsObservedIncidence(observedIncidence),
  itsDiffEvidence(itsNumAllExemplars,1),
  itsPredictedProbability(numAllExemplars(), 1),
  itsCachedLogL_1_2(0.0)
{}

Classifier::~Classifier()
{}

void Classifier::forwardProbit(double thresh, double sigmaNoise)
{
DOTRACE("Classifier::forwardProbit");

  const double divisor = (1.0 / Num::SQRT_2) * (1.0 / sigmaNoise);

  MtxConstIter diffev = itsDiffEvidence.colIter(0);

  MtxIter pp = itsPredictedProbability.colIter(0);

  for (; pp.hasMore(); ++diffev, ++pp) {
	 double alpha_val = thresh - *diffev;

	 //
	 // alpha = (thresh - diffEvidence) / sigmaNoise
	 //
	 // p = 0.5 * erfc(alpha / sqrt(2))
	 //

	 *pp = 0.5*Num::erfc(alpha_val * divisor);
  }
}

//---------------------------------------------------------------------
//
// ll = sum(gammaln(1+sum(observedIncidence,2))) - ...
//      sum(sum(gammaln(observedIncidence+1)))  + ...
//      sum(sum(observedIncidence.*log(predictedProbability)));
//
//---------------------------------------------------------------------

double Classifier::computeLogL(LogLType type)
{
DOTRACE("Classifier::computeLogL");

  if (itsCachedLogL_1_2 == 0.0)
	 {
		for(int k = 0; k < itsNumAllExemplars; ++k) {
		  double oi1 = itsObservedIncidence.at(k,0);
		  double oi2 = itsObservedIncidence.at(k,1);

		  // term 1
		  itsCachedLogL_1_2 += Num::gammaln(1.0 + oi1 + oi2);

		  // term 2
		  itsCachedLogL_1_2 -= Num::gammaln(1.0+oi1);
		  itsCachedLogL_1_2 -= Num::gammaln(1.0+oi2);
		}
	 }

  double logL_3 = 0.0;

  const double LOG_10_MINUS_50 = -115.1293;

  MtxConstIter oi1iter = itsObservedIncidence.colIter(0);
  MtxConstIter oi2iter = itsObservedIncidence.colIter(1);

  MtxConstIter ppiter = itsPredictedProbability.colIter(0);

  for(; ppiter.hasMore(); ++ppiter, ++oi1iter, ++oi2iter) {
	 double oi1 = *oi1iter;
	 double oi2 = *oi2iter;

	 // term3
	 double pp_val = (type == FULL) ? (oi1 / (oi1 + oi2)) : *ppiter;

	 if (pp_val < 1e-50) logL_3 += oi1 * LOG_10_MINUS_50;
	 else                logL_3 += oi1 * log(pp_val);

	 pp_val = 1.0 - pp_val;

	 if (pp_val < 1e-50) logL_3 += oi2 * LOG_10_MINUS_50;
	 else                logL_3 += oi2 * log(pp_val);
  }

  return itsCachedLogL_1_2 + logL_3;
}

double Classifier::currentLogL(Slice& modelParams)
{
DOTRACE("Classifier::currentLogL");

  computeDiffEv(itsObjects, modelParams, itsDiffEvidence);

  //---------------------------------------------------------------------
  //
  // Compute the predicted probabilities based on the similarity
  // differences in diffEvidence.
  //

  // the threshold is right after the DIM_OBJ_PARAMS number of
  // attentional weights
  const double thresh = modelParams[DIM_OBJ_PARAMS];

  const double sigmaNoise = computeSigmaNoise(modelParams[DIM_OBJ_PARAMS+1]);

  // predictedProbability = forwardProbit(diffEvidence, thresh, sigmaNoise);
  forwardProbit(thresh, sigmaNoise);


  //---------------------------------------------------------------------
  //
  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.
  //

  return computeLogL(CURRENT);
}

double Classifier::fullLogL()
{
DOTRACE("Classifier::fullLogL");
  return computeLogL(FULL);
}

double Classifier::deviance(Slice& modelParams)
{
DOTRACE("Classifier::deviance");

  double llc = currentLogL(modelParams);
  double llf = fullLogL();

  return -2 * (llc - llf);
}

// Count the category training exemplars
int Classifier::countCategory(int category) const
{
DOTRACE("Classifier::countCategory");
  int n = 0;
  MtxConstIter iter = itsObjCategories.colIter(0);
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

  int nobjs = countCategory(category);
  Mtx result(nobjs, DIM_OBJ_PARAMS);

  int r = 0;
  for (int i = 0; i < itsObjects.mrows(); ++i)
	 {
		if (exemplarCategory(i) == category)
		  result.row(r++) = exemplar(i);
	 }

  return result;
}

int Classifier::exemplarCategory(int i) const {
  return int(itsObjCategories.at(i));
}

Slice Classifier::exemplar(int i) const {
  return itsObjects.row(i);
}

static const char vcid_classifier_cc[] = "$Header$";
#endif // !CLASSIFIER_CC_DEFINED
