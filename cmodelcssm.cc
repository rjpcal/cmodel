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
//    itsScaledWeights1(0,0),
//    itsScaledWeights2(0,0)
{
DOTRACE("CModelCssm::CModelCssm");
}

CModelCssm::~CModelCssm() {}

void CModelCssm::loadModelParams(Slice& modelParams)
{
DOTRACE("CModelCssm::loadModelParams");

//    int nwts = numStoredExemplars() * numTrainingExemplars();

  itsScaledWeights = Mtx(modelParams);
//    itsScaledWeights1 = Mtx(modelParams.leftmost(nwts));
//    itsScaledWeights2 = Mtx(modelParams.rightmost(nwts));

  itsScaledWeights.reshape(2*numStoredExemplars(),
  									numTrainingExemplars());
//    itsScaledWeights1.reshape(numStoredExemplars(),
//  									 numTrainingExemplars());
//    itsScaledWeights2.reshape(numStoredExemplars(),
//  									 numTrainingExemplars());

  //
  // Rescale the stored exemplar weights so that they sum to 1.
  //

  itsScaledWeights.apply(abs);
//    itsScaledWeights1.apply(abs);
//    itsScaledWeights2.apply(abs);

  for (int r = 0; r < itsScaledWeights.mrows(); ++r)
//    for (int r = 0; r < itsScaledWeights1.mrows(); ++r)
	 {
		Slice row = itsScaledWeights.row(r);
		row /= row.sum();
//  		Slice row1 = itsScaledWeights1.row(r);
//  		row1 /= row1.sum();

//  		Slice row2 = itsScaledWeights2.row(r);
//  		row2 /= row2.sum();
	 }

  {for (int n = 0; n < itsStored1.mrows(); ++n)
	 {
		Slice result1(itsStored1.row(n));
  		training1().leftMultAndAssign(itsScaledWeights.row(n),
//  		training1().leftMultAndAssign(itsScaledWeights1.row(n),
												result1);
	 }
  }

  {for (int n = 0; n < itsStored2.mrows(); ++n)
	 {
		Slice result2(itsStored2.row(n));
  		training2().leftMultAndAssign(itsScaledWeights.row(n+numStoredExemplars()),
//  		training2().leftMultAndAssign(itsScaledWeights2.row(n),
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
