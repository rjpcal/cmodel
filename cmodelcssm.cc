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
  itsStored2(numStoredExemplars, DIM_OBJ_PARAMS)
{
DOTRACE("CModelCssm::CModelCssm");
}

CModelCssm::~CModelCssm() {}

void CModelCssm::loadModelParams(Slice& modelParams)
{
DOTRACE("CModelCssm::loadModelParams");

  int nex = numStoredExemplars(); 

  Mtx allScaledWeights = Mtx(modelParams);

  allScaledWeights.reshape(2*nex, numTrainingExemplars());

  Mtx scaledWeights1 = allScaledWeights.rows(0,nex);
  Mtx scaledWeights2 = allScaledWeights.rows(nex,nex);


  //
  // Rescale the stored exemplar weights so that they sum to 1.
  //

  scaledWeights1.apply(abs);
  scaledWeights2.apply(abs);

  for (int r = 0; r < scaledWeights1.mrows(); ++r)
	 {
  		Slice row1 = scaledWeights1.row(r);
  		row1 /= row1.sum();

  		Slice row2 = scaledWeights2.row(r);
  		row2 /= row2.sum();
	 }

  itsStored1.assign_MMmul(scaledWeights1, training1());
  itsStored2.assign_MMmul(scaledWeights2, training2());
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
