///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:42:01 2001
// written: Fri Mar  9 17:34:05 2001
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
  CModelExemplar(objParams, observedIncidence, countCategory(objParams,0))
{
DOTRACE("CModelGcm::CModelGcm");
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
  		return Slice(training1()[n], numAllExemplars());
	 }

  else if (CAT2 == cat)
	 {
  		return Slice(training2()[n], numAllExemplars());
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return Slice(0,0); // can't happen, but placate the compiler
}

static const char vcid_cmodelgcm_cc[] = "$Header$";
#endif // !CMODELGCM_CC_DEFINED
