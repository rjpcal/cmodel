///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:42:01 2001
// written: Fri Mar  9 17:05:24 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_CC_DEFINED
#define CMODELGCM_CC_DEFINED

#include "cmodelgcm.h"

#include "error.h"
#include "rutil.h"

#include "trace.h"

CModelGcm::CModelGcm(const Rat& objParams,
							const Rat& observedIncidence) :
  CModelExemplar(objParams, observedIncidence, countCategory(objParams,0)),
  itsNumTrainingExemplars(numStoredExemplars()),
  itsCat1(numStoredExemplars()),
  itsCat2(numStoredExemplars())
{
DOTRACE("CModelGcm::CModelGcm");
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

CModelGcm::~CModelGcm() {}

void CModelGcm::loadModelParams(Rat& modelParams)
{
DOTRACE("CModelGcm::loadModelParams");
}

Slice CModelGcm::findStoredExemplar(Category cat, int n)
{
DOTRACE("CModelGcm::findStoredExemplar");

  if (CAT1 == cat)
	 {
  		return Slice(itsCat1[n], numAllExemplars());
	 }

  else if (CAT2 == cat)
	 {
  		return Slice(itsCat2[n], numAllExemplars());
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return Slice(0,0); // can't happen, but placate the compiler
}

static const char vcid_cmodelgcm_cc[] = "$Header$";
#endif // !CMODELGCM_CC_DEFINED
