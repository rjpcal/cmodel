///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:32:31 2001
// written: Fri Mar  9 14:34:20 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_CC_DEFINED
#define CMODELEXEMPLAR_CC_DEFINED

#include "cmodelexemplar.h"

#include "rutil.h"
#include "trace.h"

#include <cmath>

double minkDist(const double* w, int nelems,
                const double* x1, int stride1,
                const double* x2, int stride2,
                double r, double r_inv)
{
  double wt_sum = 0.0;
  for (int k = 0; k < nelems; ++k)
	 {
		wt_sum +=
		  w[k] *
		  pow( abs( x1[k*stride1] - x2[k*stride2]), r);
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
  MinkDist2Binder(const double* attWeights, int nelems,
						int stride1,
						const double* x2) :
	 itsAttWeights(attWeights),
	 itsNelems(nelems),
	 itsStride1(stride1),
	 itsX2(x2)
  {}

  // Specialized Minkowski distance for r==2
  double minkDist2(const double* x1) const
  {
	 double wt_sum = 0.0;
	 const double* x2 = itsX2;
	 const double* w = itsAttWeights;
	 for (int k = 0; k < itsNelems; ++k, x1 += itsStride1, ++x2, ++w)
		{
		  wt_sum +=
			 (*w) * 
			 ((*x1) - (*x2)) * ((*x1) - (*x2));
		}
	 return sqrt(wt_sum);	 
  }

private:
  const double* const itsAttWeights;
  int itsNelems;
  int itsStride1;
  const double* const itsX2;
};


///////////////////////////////////////////////////////////////////////
//
// CModelExemplar member definitions
//
///////////////////////////////////////////////////////////////////////

CModelExemplar::CModelExemplar(const Rat& objParams,
										 const Rat& observedIncidence,
										 int numStoredExemplars) :
  Classifier(objParams,
				 observedIncidence),
  itsObjParams(objParams),
  itsNumStoredExemplars(numStoredExemplars)
{
}

CModelExemplar::~CModelExemplar()
{
}


void CModelExemplar::doDiffEvidence(const double* attWeights,
												const double* storedExemplar1,
												const double* storedExemplar2,
												double minkPower,
												double minkPowerInv)
{
DOTRACE("CModelExemplar::doDiffEvidence");

  for (int y = 0; y < numAllExemplars(); ++y) {

	 // compute similarity of ex-y to stored-1-x
	 double sim1 =
		minkDist(attWeights, DIM_OBJ_PARAMS,
					itsObjParams.data()+y+numAllExemplars(), numAllExemplars(),
					storedExemplar1, 1,
					minkPower, minkPowerInv);

	 diffEvidence(y) += exp(-sim1);

		// compute similarity of ex-y to stored-2-x
	 double sim2 =
		minkDist(attWeights, DIM_OBJ_PARAMS,
					itsObjParams.data()+y+numAllExemplars(), numAllExemplars(),
					storedExemplar2, 1,
					minkPower, minkPowerInv);

	 diffEvidence(y) -= exp(-sim2);
  }
}


void CModelExemplar::doDiffEvidence2(const double* attWeights,
												 const double* storedExemplar1,
												 const double* storedExemplar2)
{
DOTRACE("CModelExemplar::doDiffEvidence2");

  MinkDist2Binder binder1(attWeights, DIM_OBJ_PARAMS,
								  numAllExemplars(),
								  storedExemplar1);

  MinkDist2Binder binder2(attWeights, DIM_OBJ_PARAMS,
								  numAllExemplars(),
								  storedExemplar2);

  // This finds the first testExemplar data point (skipping the
  // first column of itsObjParams, which contains category labels)
  const double* testExemplar = itsObjParams.data()+numAllExemplars();

  for (int y = 0; y < numAllExemplars(); ++y, ++testExemplar) {

	 // compute similarity of ex-y to stored-1-x
	 const double sim1 = binder1.minkDist2(testExemplar);

	 // compute similarity of ex-y to stored-2-x
	 const double sim2 = binder2.minkDist2(testExemplar);

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

void CModelExemplar::computeDiffEv(Rat& modelParams) {
DOTRACE("CModelExemplar::computeDiffEv");

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  double* attWeights = modelParams.data();

  for (int i = 0; i < DIM_OBJ_PARAMS; ++i)
	 attWeights[i] = abs(attWeights[i]);

  loadModelParams(modelParams);

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  const double* stored1;
  const double* stored2;

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

	 stored1 = findStoredExemplar(CAT1, x);
	 stored2 = findStoredExemplar(CAT2, x);

	 if (minkPower == 2.0) {
		doDiffEvidence2(attWeights, stored1, stored2);
	 }
	 else {
		doDiffEvidence(attWeights, stored1, stored2,
							minkPower, minkPowerInv);
	 }
  }
}

double CModelExemplar::sigmaScalingFactor() const
{
  return sqrt(itsNumStoredExemplars*2);
}

static const char vcid_cmodelexemplar_cc[] = "$Header$";
#endif // !CMODELEXEMPLAR_CC_DEFINED
