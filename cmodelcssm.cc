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
  itsStored1(numStoredExemplars, DIM_OBJ_PARAMS),
  itsStored2(numStoredExemplars, DIM_OBJ_PARAMS),
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

  {for (int n = 0; n < itsStored1.mrows(); ++n)
	 {
		Slice result1(itsStored1.row(n));
		training1().leftMultAndAssign(itsScaledWeights.row(n),
												result1);
	 }
  }

  {for (int n = 0; n < itsStored2.mrows(); ++n)
	 {
		Slice result2(itsStored2.row(n));
		training2().leftMultAndAssign(itsScaledWeights.row(n+numStoredExemplars()),
												result2);
	 }
  }
}

const Mtx& CModelCssm::getStoredExemplars(Category cat)
{
  if (CAT1 == cat)
	 {
		return itsStored1;
	 }

  else if (CAT2 == cat)
	 {
		return itsStored2;
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return Mtx::emptyMtx(); // can't happen, but placate the compiler
}

static const char vcid_cmodelcssm_cc[] = "$Header$";
#endif // !CMODELCSSM_CC_DEFINED
