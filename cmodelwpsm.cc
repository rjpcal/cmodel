///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 17:35:56 2001
// written: Tue Mar 13 16:25:47 2001
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
  CModelExemplar(objParams, observedIncidence, 1)
{
DOTRACE("CModelWpsm::CModelWpsm");

  for (int i = 0; i < DIM_OBJ_PARAMS; ++i)
	 {
		itsPrototype1[i] = 0.0;

		int mrows1 = training1().mrows();

		{for (int k = 0; k < mrows1; ++k)
		  itsPrototype1[i] += training1().at(k,i);
		}

		itsPrototype1[i] /= mrows1;

		itsPrototype2[i] = 0.0;

		int mrows2 = training2().mrows();

		{for (int k = 0; k < mrows2; ++k)
		  itsPrototype2[i] += training2().at(k,i);
		}

		itsPrototype2[i] /= mrows2;
	 }
}


CModelWpsm::~CModelWpsm() {}


ConstSlice CModelWpsm::findStoredExemplar(Category cat, int n)
{
DOTRACE("CModelWpsm::findStoredExemplar");

  if (CAT1 == cat)
	 {
  		return Slice(&itsPrototype1[0], 1, DIM_OBJ_PARAMS);
	 }

  else if (CAT2 == cat)
	 {
  		return Slice(&itsPrototype2[0], 1, DIM_OBJ_PARAMS);
	 }

  else
	 throw ErrorWithMsg("unknown category enumerator in findStoredExemplar");

  return ConstSlice(); // can't happen, but placate the compiler
}

static const char vcid_cmodelwpsm_cc[] = "$Header$";
#endif // !CMODELWPSM_CC_DEFINED
