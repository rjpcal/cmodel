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
  itsScaledWeights(0)
{
DOTRACE("CModelCssm::CModelCssm");
}

CModelCssm::~CModelCssm() {}

void CModelCssm::scaleWeights(double* weights, int numRawWeights)
{
DOTRACE("CModelCssm::scaleWeights");

  int mrows = numStoredExemplars()*2;
  int ncols = numTrainingExemplars();

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
		(2*numTrainingExemplars()*numStoredExemplars() + 6)) {
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

Slice CModelCssm::findStoredExemplar(Category cat, int n)
{
  if (CAT1 == cat)
	 {
		Num::linearCombo(numTrainingExemplars(),
							  itsScaledWeights+n,
							  2*numStoredExemplars(),
							  &training1()[0], DIM_OBJ_PARAMS,
							  &itsStored1[0]);

		return Slice(&itsStored1[0], 1);
	 }

  else if (CAT2 == cat)
	 {
		Num::linearCombo(numTrainingExemplars(),
							  itsScaledWeights+n+numStoredExemplars(),
							  2*numStoredExemplars(),
							  &training2()[0], DIM_OBJ_PARAMS,
							  &itsStored2[0]);

		return Slice(&itsStored2[0], 1);
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return Slice(0,0); // can't happen, but placate the compiler
}

static const char vcid_cmodelcssm_cc[] = "$Header$";
#endif // !CMODELCSSM_CC_DEFINED
