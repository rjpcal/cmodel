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

double minkDist(const Slice& wts,
                MtxConstIter x1,
                MtxConstIter x2,
                double r, double r_inv)
{
  double wt_sum = 0.0;

  for (MtxConstIter wt = wts.begin();
		 wt.hasMore();
		 ++wt, ++x1, ++x2)
	 {
		wt_sum += (*wt) * pow( abs( *x1 - *x2), r);
	 }
  return pow(wt_sum, r_inv);
}

//
// This is just a functor that binds 4 of the 5 arguments to
// minkDist2, so that we don't have to keep passing 5 arguments to
// functions all the time in the time-critical inner loop.
//

class MinkDist2Binder {
public:
  MinkDist2Binder(MtxConstIter attWeights,
						MtxConstIter x2) :
	 itsAttWeights(attWeights),
	 itsX2(x2)
  {}

  // Specialized Minkowski distance for r==2
  double minkDist2(MtxConstIter x1) const
  {
	 double wt_sum = 0.0;
	 MtxConstIter wt = itsAttWeights;
	 MtxConstIter x2 = itsX2;

	 for (; wt.hasMore(); ++wt, ++x1, ++x2)
		{
		  wt_sum +=
			 (*wt) * 
			 ((*x1) - (*x2)) * ((*x1) - (*x2));
		}
	 return sqrt(wt_sum);	 
  }

private:
  const MtxConstIter itsAttWeights;
  const MtxConstIter itsX2;
};


///////////////////////////////////////////////////////////////////////
//
// CModelExemplar member definitions
//
///////////////////////////////////////////////////////////////////////

CModelExemplar::CModelExemplar(const Mtx& objParams,
										 const Mtx& observedIncidence,
										 int numStoredExemplars) :
  Classifier(objParams,
				 observedIncidence),
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
			 itsTraining1.row(c1++) = objParams.row(i).rightmost(DIM_OBJ_PARAMS);
		  }
		else if (int(objParams.at(i,0)) == 1)
		  {
			 itsTraining2.row(c2++) = objParams.row(i).rightmost(DIM_OBJ_PARAMS);
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

void CModelExemplar::doDiffEvidence(const Slice& attWeights,
												const Slice& storedExemplar1,
												const Slice& storedExemplar2,
												double minkPower,
												double minkPowerInv)
{
DOTRACE("CModelExemplar::doDiffEvidence");

  throw ErrorWithMsg("doDiffEvidence not implemented");

  for (int y = 0; y < numAllExemplars(); ++y) {

	 // compute similarity of ex-y to stored-1-x
	 double sim1 =
		minkDist(attWeights,
					exemplar(y).begin(),
					storedExemplar1.begin(),
					minkPower, minkPowerInv);

	 diffEvidence(y) += exp(-sim1);

	 // compute similarity of ex-y to stored-2-x
	 double sim2 =
		minkDist(attWeights,
					exemplar(y).begin(),
					storedExemplar2.begin(),
					minkPower, minkPowerInv);

	 diffEvidence(y) -= exp(-sim2);
  }
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

  if (minkPower == 2.0) {
	 MtxConstIter attWts = attWeights.begin();

	 minivec<MtxConstIter> exemplars;

	 for (int yy = 0; yy < numAllExemplars(); ++yy) {
		exemplars.push_back(exemplar(yy).begin());
	 }

	 for (int x = 0; x < itsNumStoredExemplars; ++x) {

		MinkDist2Binder binder1(attWts, stored1.rowIter(x));
		MinkDist2Binder binder2(attWts, stored2.rowIter(x));

		for (int y = 0; y < numAllExemplars(); ++y) {

		  // compute similarity of ex-y to stored-1-x
		  const double sim1 = binder1.minkDist2(exemplars[y]);

		  // compute similarity of ex-y to stored-2-x
		  const double sim2 = binder2.minkDist2(exemplars[y]);

		  diffEvidence(y) += exp(-sim1) - exp(-sim2);
		}
	 }
  }
  else {
	 for (int x = 0; x < itsNumStoredExemplars; ++x) {
		doDiffEvidence(attWeights, stored1.row(x), stored2.row(x),
							minkPower, minkPowerInv);
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
