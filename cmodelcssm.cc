///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:25:38 2001
// written: Tue Mar 13 17:58:58 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_CC_DEFINED
#define CMODELCSSM_CC_DEFINED

#include "cmodelcssm.h"

#include "error.h"
#include "num.h"
#include "mtx.h"
#include "trace.h"

#include <cmath>

CModelCssm::CModelCssm(const Mtx& objParams,
							  const Mtx& observedIncidence,
							  int numStoredExemplars) :
  CModelExemplar(objParams, observedIncidence, numStoredExemplars),
  itsStored1(1, DIM_OBJ_PARAMS),
  itsStored2(1, DIM_OBJ_PARAMS),
  itsScaledWeights(0,0)
{
DOTRACE("CModelCssm::CModelCssm");
}

CModelCssm::~CModelCssm() {}

void CModelCssm::loadModelParams(Slice& modelParams)
{
DOTRACE("CModelCssm::loadModelParams");

  itsScaledWeights = Mtx(modelParams);

  itsScaledWeights.reshape(2*numStoredExemplars(),
									numTrainingExemplars());

  //
  // Rescale the stored exemplar weights so that they sum to 1.
  //

  itsScaledWeights.apply(abs);

  for (int r = 0; r < itsScaledWeights.mrows(); ++r)
	 {
		Slice row = itsScaledWeights.row(r);
		row /= row.sum();
	 }
}

ConstSlice CModelCssm::findStoredExemplar(Category cat, int n)
{
  if (CAT1 == cat)
	 {
		Slice result(itsStored1.row(0));

		training1().leftMultAndAssign(itsScaledWeights.row(n),
  												result);

		return result;
	 }

  else if (CAT2 == cat)
	 {
  		Slice result(itsStored2.row(0));

		training2().leftMultAndAssign(itsScaledWeights.row(n+numStoredExemplars()),
  												result);

  		return result;
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return ConstSlice(); // can't happen, but placate the compiler
}

static const char vcid_cmodelcssm_cc[] = "$Header$";
#endif // !CMODELCSSM_CC_DEFINED
