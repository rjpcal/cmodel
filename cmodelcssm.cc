///////////////////////////////////////////////////////////////////////
//
// modelcssm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:25:38 2001
// written: Thu Mar  8 16:30:34 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MODELCSSM_CC_DEFINED
#define MODELCSSM_CC_DEFINED

#include "modelcssm.h"

#include "error.h"
#include "num.h"
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
// ModelCssm member definitions
//
///////////////////////////////////////////////////////////////////////

ModelCssm::ModelCssm(const Rat& objParams,
							const Rat& observedIncidence,
							int numStoredExemplars) :
  Classifier(objParams,
				 observedIncidence),
  itsObjParams(objParams),
  itsNumStoredExemplars(numStoredExemplars),
  itsNum1(countCategory(objParams, 0)),
  itsNum2(countCategory(objParams, 1)),
  itsNumTrainingExemplars(itsNum1),
  itsCat1(new constDblPtr[itsNum1]),
  itsCat2(new constDblPtr[itsNum2])
{
  if (itsNum1 != itsNum2) {
	 throw ErrorWithMsg("the two categories must have the "
							  "same number of training exemplars");
  }

  // Find the category 1 and category 2 training exemplars
  int c1=0,c2=0;
  for (int i = 0; i < itsObjParams.mrows(); ++i)
	 {
		if (int(itsObjParams.at(i,0)) == 0)
		  itsCat1[c1++] = itsObjParams.address(i,1);
		else if (int(itsObjParams.at(i,0)) == 1)
		  itsCat2[c2++] = itsObjParams.address(i,1);
	 }
}

ModelCssm::~ModelCssm()
{
  delete [] itsCat2;
  delete [] itsCat1;
}

// Count the category training exemplars
int ModelCssm::countCategory(const Rat& params, int category) {
  int n = 0;
  for (int i = 0; i < params.mrows(); ++i)
	 {
		if (int(params.at(i,0)) == category)
		  ++n;
	 }
  return n;
}


void ModelCssm::scaleWeights(double* weights, int numRawWeights)
{
DOTRACE("ModelCssm::scaleWeights");

  int mrows = itsNumStoredExemplars*2;
  int ncols = itsNumTrainingExemplars;

  if ( numRawWeights != (mrows*ncols) )
	 throw ErrorWithMsg("weights must have "
							  "2*numStoredExemplars*numTrainingExemplars elements");

  for (int i = 0; i < mrows; ++i)
	 {
		double sum_wt = 0.0;
		{
		  for (int ix = i; ix < numRawWeights; ix+=mrows)
			 sum_wt += abs(weights[ix]);
		}
		{
		  for (int ix = i; ix < numRawWeights; ix+=mrows)
			 weights[ix] = abs(weights[ix]) / sum_wt;
		}
	 }
}


void ModelCssm::computeSimilarity(const double* attWeights,
											 const double* storedExemplar1,
											 const double* storedExemplar2,
											 double minkPower,
											 double minkPowerInv)
{
DOTRACE("ModelCssm::computeSimilarity");

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


void ModelCssm::computeSimilarity2(const double* attWeights,
											  const double* storedExemplar1,
											  const double* storedExemplar2)
{
DOTRACE("ModelCssm::computeSimilarity2");

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

void ModelCssm::computeDiffEv(Rat& modelParams) {
DOTRACE("ModelCssm::computeDiffEv");

  if (modelParams.length() !=
		(2*itsNumTrainingExemplars*itsNumStoredExemplars + 6)) {
    throw ErrorWithMsg("wrong number of model parameters");
  }

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  double* attWeights = modelParams.data();

  {
    for (int i = 0; i < DIM_OBJ_PARAMS; ++i)
      attWeights[i] = abs(attWeights[i]);
  }

  //---------------------------------------------------------------------
  //
  // Rescale the stored exemplar weights so that they sum to 1.
  //

  double* rawWeights = modelParams.data() + 6;
  const int numRawWeights = modelParams.nelems() - 6;

  scaleWeights(rawWeights, numRawWeights);

  const double* scaledWeights = rawWeights;

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  double stored1[DIM_OBJ_PARAMS];
  double stored2[DIM_OBJ_PARAMS];

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

    Num::linearCombo(itsNumTrainingExemplars, scaledWeights+x,
							2*itsNumStoredExemplars,
							itsCat1, numAllExemplars(), DIM_OBJ_PARAMS,
							&stored1[0]);

    Num::linearCombo(itsNumTrainingExemplars, scaledWeights+x+itsNumStoredExemplars,
							2*itsNumStoredExemplars,
							itsCat2, numAllExemplars(), DIM_OBJ_PARAMS,
							&stored2[0]);

	 if (minkPower == 2.0) {
		computeSimilarity2(attWeights, &stored1[0], &stored2[0]);
	 }
	 else {
		computeSimilarity(attWeights, &stored1[0], &stored2[0],
								minkPower, minkPowerInv);
	 }
  }
}

double ModelCssm::sigmaScalingFactor() const
{
  return sqrt(itsNumStoredExemplars*2);
}

static const char vcid_modelcssm_cc[] = "$Header$";
#endif // !MODELCSSM_CC_DEFINED
