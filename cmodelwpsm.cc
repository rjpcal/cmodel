///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 17:35:56 2001
// written: Tue Mar 13 18:02:18 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELWPSM_CC_DEFINED
#define CMODELWPSM_CC_DEFINED

#include "cmodelwpsm.h"

#include "error.h"
#include "mtx.h"

#include "trace.h"


CModelWpsm::CModelWpsm(const Mtx& objParams,
							  const Mtx& observedIncidence) :
  CModelExemplar(objParams, observedIncidence, 1),
  itsPrototype1(1, DIM_OBJ_PARAMS),
  itsPrototype2(1, DIM_OBJ_PARAMS)
{
DOTRACE("CModelWpsm::CModelWpsm");

  for (int i = 0; i < DIM_OBJ_PARAMS; ++i)
	 {
		itsPrototype1.at(0,i) = training1().column(i).mean();
		itsPrototype2.at(0,i) = training2().column(i).mean();
	 }
}


CModelWpsm::~CModelWpsm() {}


ConstSlice CModelWpsm::findStoredExemplar(Category cat, int n)
{
DOTRACE("CModelWpsm::findStoredExemplar");

  if (CAT1 == cat)
	 {
  		return itsPrototype1.row(0);
	 }

  else if (CAT2 == cat)
	 {
  		return itsPrototype2.row(0);
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return ConstSlice(); // can't happen, but placate the compiler
}

static const char vcid_cmodelwpsm_cc[] = "$Header$";
#endif // !CMODELWPSM_CC_DEFINED
