///////////////////////////////////////////////////////////////////////
//
// classifier.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:34:12 2001
// written: Thu Mar  8 17:03:16 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_CC_DEFINED
#define CLASSIFIER_CC_DEFINED

#include "classifier.h"

#include "error.h"
#include "modelcssm.h"
#include "num.h"
#include "rutil.h"
#include "strings.h"
#include "trace.h"

#include "util/pointers.h"

///////////////////////////////////////////////////////////////////////
//
// Classifier member definitions
//
///////////////////////////////////////////////////////////////////////

Classifier::Classifier(const Rat& objParams,
							  const Rat& observedIncidence) :
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

  const double* observedIncidence1 =
	 itsObservedIncidence.data();

  const double* observedIncidence2 =
	 itsObservedIncidence.data()+itsNumAllExemplars;

  for(int k = 0; k < itsNumAllExemplars; ++k) {
	 double oi1 = observedIncidence1[k];
	 double oi2 = observedIncidence2[k];

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

shared_ptr<Classifier> Classifier::make(const fixed_string& whichType,
													 const Rat& objParams,
													 const Rat& observedIncidence,
													 int numStoredExemplars)
{
DOTRACE("Classifier::make");
  if (whichType == "cssm")
	 return shared_ptr<Classifier>(
		new ModelCssm(objParams, observedIncidence,numStoredExemplars));
  else
	 {
		ErrorWithMsg err("unknown classifier type: ");
		err.appendMsg(whichType.c_str());
		throw err;
	 }
}

double Classifier::currentLogL(Rat& modelParams)
{
DOTRACE("Classifier::currentLogL");

  resetDiffEv();

  computeDiffEv(modelParams);


  //---------------------------------------------------------------------
  //
  // Compute the predicted probabilities based on the similarity
  // differences in diffEvidence.
  //

  // thresh = modelParams(5);
  const double thresh = modelParams.at(4);

  // sigmaNoise = sqrt(n1+n2)*modelParams(6);
  const double sigmaNoise = sigmaScalingFactor() * modelParams.at(5);

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

double Classifier::deviance(Rat& modelParams) {
  double llc = currentLogL(modelParams);
  double llf = fullLogL();

  return -2 * (llc - llf);
}

static const char vcid_classifier_cc[] = "$Header$";
#endif // !CLASSIFIER_CC_DEFINED
