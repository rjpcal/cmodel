///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:34:12 2001
// written: Mon Mar 12 14:53:29 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

#include "classifier.h"

#include "error.h"
#include "num.h"
#include "mtx.h"
#include "strings.h"
#include "trace.h"

#include "util/pointers.h"

///////////////////////////////////////////////////////////////////////
//
// Classifier member definitions
//
///////////////////////////////////////////////////////////////////////

Classifier::Classifier(const Mtx& objParams,
							  const Mtx& observedIncidence) :
  itsObjParams(objParams),
  itsNumAllExemplars(objParams.mrows()),
  itsObservedIncidence(observedIncidence),
  itsDiffEvidence(new double[itsNumAllExemplars]),
  itsPredictedProbability(new double[numAllExemplars()])
{}

Classifier::~Classifier()
{
  delete [] itsPredictedProbability;
  delete [] itsDiffEvidence;
}

void Classifier::forwardProbit(double thresh, double sigmaNoise) const
{
DOTRACE("Classifier::forwardProbit");

  const double divisor = (1.0 / Num::SQRT_2) * (1.0 / sigmaNoise);

  for (int i = 0; i < itsNumAllExemplars; ++i) {
	 double alpha_val = thresh-itsDiffEvidence[i];

	 //
	 // alpha = (thresh - diffEvidence) / sigmaNoise
	 //
	 // p = 0.5 * erfc(alpha / sqrt(2))
	 //

	 itsPredictedProbability[i] = 0.5*Num::erfc(alpha_val * divisor);
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
  double ll = 0.0;

  const double LOG_10_MINUS_50 = -115.1293;

  for(int k = 0; k < itsNumAllExemplars; ++k) {
	 double oi1 = itsObservedIncidence.at(k,0);
	 double oi2 = itsObservedIncidence.at(k,1);

	 // term 1
	 ll += Num::gammaln(1.0 + oi1 + oi2);

	 // term 2
	 ll -= Num::gammaln(1.0+oi1);
	 ll -= Num::gammaln(1.0+oi2);

	 // term3
	 double pp_val = (type == FULL) ?
		(oi1 / (oi1 + oi2)) :
		itsPredictedProbability[k];

	 if (pp_val < 1e-50) ll += oi1 * LOG_10_MINUS_50;
	 else                ll += oi1 * log(pp_val);

	 pp_val = 1.0 - pp_val;

	 if (pp_val < 1e-50) ll += oi2 * LOG_10_MINUS_50;
	 else                ll += oi2 * log(pp_val);
  }

  return ll;
}

double Classifier::currentLogL(Mtx& modelParams)
{
DOTRACE("Classifier::currentLogL");

  resetDiffEv();

  computeDiffEv(modelParams);


  //---------------------------------------------------------------------
  //
  // Compute the predicted probabilities based on the similarity
  // differences in diffEvidence.
  //

  // the threshold is right after the DIM_OBJ_PARAMS number of
  // attentional weights
  const double thresh = modelParams.at(DIM_OBJ_PARAMS);

  const double sigmaNoise = fetchSigmaNoise(modelParams);

  // predictedProbability = forwardProbit(diffEvidence, thresh, sigmaNoise);
  forwardProbit(thresh, sigmaNoise);


  //---------------------------------------------------------------------
  //
  // Compute the loglikelihood based on the predicted probabilities
  // and the observed incidences.
  //

  return computeLogL(CURRENT);
}

double Classifier::fullLogL() {
  return computeLogL(FULL);
}

double Classifier::deviance(Mtx& modelParams) {
  double llc = currentLogL(modelParams);
  double llf = fullLogL();

  return -2 * (llc - llf);
}

int Classifier::exemplarCategory(int i) const {
  return int(itsObjParams.at(i));
}

Slice Classifier::exemplar(int i) const {
  // Skip the first column which contains category info
  return itsObjParams.rowSlice(i).rightmost(DIM_OBJ_PARAMS);
}

static const char vcid_classifier_cc[] = "$Header$";
#endif // !CLASSIFIER_CC_DEFINED
