///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:32:31 2001
// written: Mon Mar 12 17:05:30 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_CC_DEFINED
#define CMODELEXEMPLAR_CC_DEFINED

#include "cmodelexemplar.h"

#include "error.h"
#include "mtx.h"
#include "trace.h"

#include <cmath>

double minkDist(const ConstSlice& wts,
                ConstSlice::ConstIterator x1,
                ConstSlice::ConstIterator x2,
                double r, double r_inv)
{
  double wt_sum = 0.0;

  for (ConstSlice::ConstIterator
			wt = wts.begin(),
			wend = wts.end();
		 wt != wend;
		 ++wt)
	 {
		wt_sum += (*wt) * pow( abs( *x1 - *x2), r);

		++wt;
		++x1;
		++x2;
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
  MinkDist2Binder(const ConstSlice& attWeights,
						const ConstSlice& x2) :
	 itsAttWeights(attWeights),
	 itsNelems(attWeights.nelems()),
	 itsX2(x2)
  {}

  // Specialized Minkowski distance for r==2
  double minkDist2(Slice::ConstIterator x1) const
  {
	 double wt_sum = 0.0;
	 Slice::ConstIterator wt = itsAttWeights.begin();
	 Slice::ConstIterator x2 = itsX2.begin();

	 for (int k = 0; k < itsNelems; ++k, ++x1, ++x2, ++wt)
		{
		  wt_sum +=
			 (*wt) * 
			 ((*x1) - (*x2)) * ((*x1) - (*x2));
		}
	 return sqrt(wt_sum);	 
  }

private:
  const ConstSlice itsAttWeights;
  int itsNelems;
  const ConstSlice itsX2;
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
  itsTraining1(itsNumTrainingExemplars),
  itsTraining2(itsNumTrainingExemplars),
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
		  itsTraining1[c1++].rebind(
						objParams.row(i).rightmost(DIM_OBJ_PARAMS));
		else if (int(objParams.at(i,0)) == 1)
		  itsTraining2[c2++].rebind(
						objParams.row(i).rightmost(DIM_OBJ_PARAMS));
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

void CModelExemplar::doDiffEvidence(const ConstSlice& attWeights,
												const ConstSlice& storedExemplar1,
												const ConstSlice& storedExemplar2,
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


void CModelExemplar::doDiffEvidence2(const ConstSlice& attWeights,
												 const ConstSlice& storedExemplar1,
												 const ConstSlice& storedExemplar2)
{
DOTRACE("CModelExemplar::doDiffEvidence2");

  MinkDist2Binder binder1(attWeights,
								  storedExemplar1);

  MinkDist2Binder binder2(attWeights,
								  storedExemplar2);

  for (int y = 0; y < numAllExemplars(); ++y) {

	 ConstSlice ex(exemplar(y));

	 // compute similarity of ex-y to stored-1-x
	 const double sim1 = binder1.minkDist2(ex.begin());

	 // compute similarity of ex-y to stored-2-x
	 const double sim2 = binder2.minkDist2(ex.begin());

	 diffEvidence(y) += exp(-sim1) - exp(-sim2);
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

void CModelExemplar::computeDiffEv(Mtx& modelParams) {
DOTRACE("CModelExemplar::computeDiffEv");

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  Slice attWeights = modelParams.asSlice().leftmost(DIM_OBJ_PARAMS);

  for (int i = 0; i < attWeights.nelems(); ++i)
	 attWeights[i] = abs(attWeights[i]);

  Slice otherParams = modelParams.asSlice().rightmost(modelParams.nelems()-6);

  loadModelParams(otherParams);

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

	 ConstSlice stored1(findStoredExemplar(CAT1, x));
	 ConstSlice stored2(findStoredExemplar(CAT2, x));

	 if (minkPower == 2.0) {
		doDiffEvidence2(attWeights, stored1, stored2);
	 }
	 else {
		doDiffEvidence(attWeights, stored1, stored2,
							minkPower, minkPowerInv);
	 }
  }
}

double CModelExemplar::fetchSigmaNoise(const Mtx& modelParams) const
{
  return modelParams.at(5) * sqrt(itsNumStoredExemplars*2.0);
}

void CModelExemplar::loadModelParams(Slice& modelParams) {}

static const char vcid_cmodelexemplar_cc[] = "$Header$";
#endif // !CMODELEXEMPLAR_CC_DEFINED
