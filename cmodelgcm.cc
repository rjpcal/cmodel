///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:42:01 2001
// written: Mon Mar 12 17:06:05 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_CC_DEFINED
#define CMODELGCM_CC_DEFINED

#include "cmodelgcm.h"

#include "error.h"
#include "mtx.h"

#include "trace.h"


CModelGcm::CModelGcm(const Mtx& objParams,
							const Mtx& observedIncidence) :
  CModelExemplar(objParams, observedIncidence, countCategory(objParams,0))
{
DOTRACE("CModelGcm::CModelGcm");
}


CModelGcm::~CModelGcm() {}


ConstSlice CModelGcm::findStoredExemplar(Category cat, int n)
{
DOTRACE("CModelGcm::findStoredExemplar");

  if (CAT1 == cat)
	 {
  		return training1()[n];
	 }

  else if (CAT2 == cat)
	 {
  		return training2()[n];
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return Slice(); // can't happen, but placate the compiler
}

static const char vcid_cmodelgcm_cc[] = "$Header$";
#endif // !CMODELGCM_CC_DEFINED
