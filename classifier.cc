///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:34:12 2001
// written: Wed Apr 18 15:43:44 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

#include "classifier.h"

#include "error.h"
#include "mxwrapper.h"
#include "num.h"
#include "strings.h"

#include "trace.h"
#include "debug.h"

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
  itsCachedLogL_1_2(0.0)
{}

Classifier::~Classifier()
{}

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

double Classifier::computeLogL(const Mtx& predictedProbability)
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

  MtxConstIter oi1iter = itsObservedIncidence.columnIter(0);
  MtxConstIter oi2iter = itsObservedIncidence.columnIter(1);

  MtxConstIter ppiter = predictedProbability.columnIter(0);

  for(; ppiter.hasMore(); ++ppiter, ++oi1iter, ++oi2iter) {
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
Mtx Classifier::classifyObjects(Slice& modelParams, const Mtx& objects)
{
DOTRACE("Classifier::classifyObjects");

  itsDiffEvidence.resize(objects.mrows(), 1);

  computeDiffEv(objects, modelParams, itsDiffEvidence);

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

double Classifier::currentLogL(Slice& modelParams)
{
DOTRACE("Classifier::currentLogL");

  Mtx pp = classifyObjects(modelParams, itsObjects);

  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.

  return computeLogL(pp);
}

double Classifier::fullLogL()
{
DOTRACE("Classifier::fullLogL");

  Mtx observedProb(numAllExemplars(), 1);
  MtxIter opiter = observedProb.columnIter(0);

  MtxConstIter oi1iter = itsObservedIncidence.columnIter(0);
  MtxConstIter oi2iter = itsObservedIncidence.columnIter(1);

  for (; opiter.hasMore(); ++opiter, ++oi1iter, ++oi2iter)
	 *opiter = (*oi1iter / (*oi1iter + *oi2iter));

  return computeLogL(observedProb);
}

double Classifier::deviance(Slice& modelParams)
{
DOTRACE("Classifier::deviance");

  double llc = currentLogL(modelParams);
  double llf = fullLogL();

  return -2 * (llc - llf);
}

Classifier::RequestResult Classifier::handleRequest(fixed_string action,
																	 const Mtx& allModelParams,
																	 mxArray* extraArgs_mx)
{
DOTRACE("Classifier::handleRequest");

  int multiplier = 1;
  // check for minus sign
  if (action.c_str()[0] == '-')
	 {
		multiplier = -1;
		fixed_string trimmed = action.c_str() + 1;
		action = trimmed;
	 }


  //---------------------------------------------------------------------
  //
  // Call the requested computational function for each set of model
  // params
  //
  //---------------------------------------------------------------------

  if ( action == "ll" || action == "llc" )
	 {
		Mtx result(allModelParams.ncols(), 1);
		for (int i = 0; i < allModelParams.ncols(); ++i)
		  {
			 Slice modelParams(allModelParams.column(i));
			 result.at(i) = multiplier * this->currentLogL(modelParams);
		  }

		return result;
	 }

  if ( action == "llf" )
	 {
		Mtx result(allModelParams.ncols(), 1);
		for (int i = 0; i < allModelParams.ncols(); ++i)
		  {
			 result.at(i) = multiplier * this->fullLogL();
		  }

		return result;
	 }

  if ( action == "dev" )
	 {
		Mtx result(allModelParams.ncols(), 1);
		for (int i = 0; i < allModelParams.ncols(); ++i)
		  {
			 Slice modelParams(allModelParams.column(i));
			 result.at(i) = multiplier * this->deviance(modelParams);
		  }

		return result;
	 }

  if ( action == "classify" )
	 {
		Mtx testObjects =
		  Mtx(MxWrapper::extractStructField(extraArgs_mx, "testObjects"),
				Mtx::COPY);

		Mtx result(testObjects.mrows(), allModelParams.ncols());

		for (int i = 0; i < allModelParams.ncols(); ++i)
		  {
			 Slice modelParams(allModelParams.column(i));
			 result.column(i) = this->classifyObjects(modelParams, testObjects);
		  }

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
