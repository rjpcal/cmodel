///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:25:38 2001
// written: Fri Apr  6 10:27:20 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_CC_DEFINED
#define CMODELCSSM_CC_DEFINED

#include "cmodelcssm.h"

#include "error.h"
#include "num.h"
#include "mtx.h"

#include <cmath>

#include "trace.h"

CModelCssm::CModelCssm(const Mtx& objParams,
							  const Mtx& observedIncidence,
							  TransferFunction transferFunc,
							  int numStoredExemplars) :
  CModelExemplar(objParams, observedIncidence,
					  numStoredExemplars, transferFunc),
  itsStored1(numStoredExemplars, DIM_OBJ_PARAMS),
  itsStored2(numStoredExemplars, DIM_OBJ_PARAMS),
  itsCachedRawWts1(numStoredExemplars, numTrainingExemplars()),
  itsCachedRawWts2(numStoredExemplars, numTrainingExemplars())
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

  bool do1 = (scaledWeights1 != itsCachedRawWts1);
  bool do2 = (scaledWeights2 != itsCachedRawWts2);

  if (do1) itsCachedRawWts1 = scaledWeights1; itsCachedRawWts1.makeUnique();
  if (do2) itsCachedRawWts2 = scaledWeights2; itsCachedRawWts2.makeUnique();

  if (do1) scaledWeights1.apply(abs);
  if (do2) scaledWeights2.apply(abs);

  for (int r = 0; r < scaledWeights1.mrows(); ++r)
	 {
		if (do1) { Slice row1 = scaledWeights1.row(r); row1 /= row1.sum(); }
		if (do2) { Slice row2 = scaledWeights2.row(r); row2 /= row2.sum(); }
	 }

  if (do1) itsStored1.assign_MMmul(scaledWeights1, training1());
  if (do2) itsStored2.assign_MMmul(scaledWeights2, training2());
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
