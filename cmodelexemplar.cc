///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:32:31 2001
// written: Thu Mar 15 15:51:15 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_CC_DEFINED
#define CMODELEXEMPLAR_CC_DEFINED

#include "cmodelexemplar.h"

#include "error.h"
#include "mtx.h"
#include "trace.h"

#include "minivec.h"

#include <cmath>

//
// This is just a functor that binds arguments to minkDist, so that
// we don't have to keep passing extra arguments to functions all the
// time in the time-critical inner loop.
//

class MinkDistBinder {
public:
  MinkDistBinder(MtxConstIter attWeights, MtxConstIter x2,
					  double r = 2.0, double r_inv = 0.5) :
	 itsAttWeights(attWeights),
	 itsX2(x2),
	 itsR(r),
	 itsRinv(r_inv)
  {}

  double minkDist(MtxConstIter x1) const
  {
	 double wt_sum = 0.0;
	 MtxConstIter wt = itsAttWeights;
	 MtxConstIter x2 = itsX2;

	 for (; wt.hasMore(); ++wt, ++x1, ++x2)
		{
		  wt_sum += (*wt) * pow( abs( *x1 - *x2), itsR);
		}
	 return pow(wt_sum, itsRinv);
  }

  // Specialized Minkowski distance for r==2
  double minkDist2(MtxConstIter x1) const
  {
	 double wt_sum = 0.0;
	 MtxConstIter wt = itsAttWeights;
	 MtxConstIter x2 = itsX2;

	 for (; wt.hasMore(); ++wt, ++x1, ++x2)
		{
		  wt_sum += (*wt) * ((*x1) - (*x2)) * ((*x1) - (*x2));
		}
	 return sqrt(wt_sum);	 
  }

private:
  const MtxConstIter itsAttWeights;
  const MtxConstIter itsX2;
  const double itsR;
  const double itsRinv;
};


///////////////////////////////////////////////////////////////////////
//
// CModelExemplar member definitions
//
///////////////////////////////////////////////////////////////////////

CModelExemplar::CModelExemplar(const Mtx& objParams,
										 const Mtx& observedIncidence,
										 int numStoredExemplars) :
  Classifier(objParams, observedIncidence),
  itsNumTrainingExemplars(countCategory(objParams,0)),
  itsTraining1(itsNumTrainingExemplars, DIM_OBJ_PARAMS),
  itsTraining2(itsNumTrainingExemplars, DIM_OBJ_PARAMS),
  itsNumStoredExemplars(numStoredExemplars)
{
  int num2 = countCategory(objParams, 1);

  if (itsNumTrainingExemplars != num2) {
	 throw ErrorWithMsg("the two categories must have the "
							  "same number of training exemplars");
  }

  // Find the category 1 and category 2 training exemplars
  int c1=0,c2=0;
  for (int i = 0; i < objParams.mrows(); ++i)
	 {
		if (int(objParams.at(i,0)) == 0)
		  {
			 itsTraining1.row(c1++) = exemplar(i);
		  }
		else if (int(objParams.at(i,0)) == 1)
		  {
			 itsTraining2.row(c2++) = exemplar(i);
		  }
	 }
}

CModelExemplar::~CModelExemplar()
{}

// Count the category training exemplars
int CModelExemplar::countCategory(const Mtx& params, int category) {
  int n = 0;
  for (int i = 0; i < params.mrows(); ++i)
	 {
		if (int(params.at(i,0)) == category)
		  ++n;
	 }
  return n;
}

//---------------------------------------------------------------------
//
// compute the minus loglikelihood for the constrained summed
// similarity model (cssm)
//
// based on the MATLAB function llcssm 
//
//     function ll = llcssm(modelParams, objParams, ...
//                          observedIncidence, numStoredExemplars)
//
//---------------------------------------------------------------------

void CModelExemplar::computeDiffEv(Slice& modelParams) {
DOTRACE("CModelExemplar::computeDiffEv");

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  Slice attWeights = modelParams.leftmost(DIM_OBJ_PARAMS);

  attWeights.apply(abs);

  Slice otherParams = modelParams.rightmost(modelParams.nelems()-
														  (DIM_OBJ_PARAMS+2));

  loadModelParams(otherParams);

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  const Mtx& stored1 = getStoredExemplars(CAT1);
  const Mtx& stored2 = getStoredExemplars(CAT2);

  MtxConstIter attWts = attWeights.begin();

  minivec<MtxConstIter> exemplars;

  for (int yy = 0; yy < numAllExemplars(); ++yy) {
	 exemplars.push_back(exemplar(yy).begin());
  }

  double sim1, sim2;

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

	 MinkDistBinder binder1(attWts, stored1.rowIter(x),
									minkPower, minkPowerInv);
	 MinkDistBinder binder2(attWts, stored2.rowIter(x),
									minkPower, minkPowerInv);

	 for (int y = 0; y < numAllExemplars(); ++y) {

		if (minkPower == 2.0) {
		  sim1 = binder1.minkDist2(exemplars[y]);
		  sim2 = binder2.minkDist2(exemplars[y]);
		}
		else {
		  sim1 = binder1.minkDist(exemplars[y]);
		  sim2 = binder2.minkDist(exemplars[y]);
		}
		diffEvidence(y) += exp(-sim1) - exp(-sim2);
	 }
  }
}

double CModelExemplar::computeSigmaNoise(double rawSigma) const
{
  return rawSigma * sqrt(itsNumStoredExemplars*2.0);
}

void CModelExemplar::loadModelParams(Slice& modelParams) {}

static const char vcid_cmodelexemplar_cc[] = "$Header$";
#endif // !CMODELEXEMPLAR_CC_DEFINED
