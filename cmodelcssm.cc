///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:25:38 2001
// written: Fri Mar  9 14:34:46 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_CC_DEFINED
#define CMODELCSSM_CC_DEFINED

#include "cmodelcssm.h"

#include "error.h"
#include "num.h"
#include "rutil.h"
#include "trace.h"

#include <cmath>

CModelCssm::CModelCssm(const Rat& objParams,
							  const Rat& observedIncidence,
							  int numStoredExemplars) :
  CModelExemplar(objParams, observedIncidence, numStoredExemplars),
  itsNumTrainingExemplars(countCategory(objParams,0)),
  itsCat1(itsNumTrainingExemplars),
  itsCat2(itsNumTrainingExemplars),
  itsScaledWeights(0)
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
		  itsCat1[c1++] = objParams.address(i,1);
		else if (int(objParams.at(i,0)) == 1)
		  itsCat2[c2++] = objParams.address(i,1);
	 }
}

CModelCssm::~CModelCssm() {}

void CModelCssm::scaleWeights(double* weights, int numRawWeights)
{
DOTRACE("CModelCssm::scaleWeights");

  int mrows = numStoredExemplars()*2;
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

void CModelCssm::loadModelParams(Rat& modelParams)
{
DOTRACE("CModelCssm::loadModelParams");

  if (modelParams.length() !=
		(2*itsNumTrainingExemplars*numStoredExemplars() + 6)) {
    throw ErrorWithMsg("wrong number of model parameters");
  }

  //
  // Rescale the stored exemplar weights so that they sum to 1.
  //

  double* rawWeights = modelParams.data() + 6;
  const int numRawWeights = modelParams.nelems() - 6;

  scaleWeights(rawWeights, numRawWeights);

  itsScaledWeights = rawWeights;
}

const double* CModelCssm::findStoredExemplar(Category cat, int n)
{
  if (CAT1 == cat)
	 {
		Num::linearCombo(itsNumTrainingExemplars,
							  itsScaledWeights+n,
							  2*numStoredExemplars(),
							  &itsCat1[0], numAllExemplars(), DIM_OBJ_PARAMS,
							  &itsStored1[0]);

		return &itsStored1[0];
	 }

  else if (CAT2 == cat)
	 {
		Num::linearCombo(itsNumTrainingExemplars,
							  itsScaledWeights+n+numStoredExemplars(),
							  2*numStoredExemplars(),
							  &itsCat2[0], numAllExemplars(), DIM_OBJ_PARAMS,
							  &itsStored2[0]);

		return &itsStored2[0];
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return 0; // can't happen, but placate the compiler
}

static const char vcid_cmodelcssm_cc[] = "$Header$";
#endif // !CMODELCSSM_CC_DEFINED
