///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 17:35:56 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELWPSM_CC_DEFINED
#define CMODELWPSM_CC_DEFINED

#include "cmodel/cmodelwpsm.h"

#include "mtx/mtx.h"

#include "util/error.h"

#include "util/trace.h"


CModelWpsm::CModelWpsm(const Mtx& objParams,
                       TransferFunction transferFunc) :
  CModelExemplar(objParams, 1, transferFunc),
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


const Mtx& CModelWpsm::getStoredExemplars(Category cat)
{
DOTRACE("CModelWpsm::findStoredExemplar");

  if (CAT1 == cat)
    {
      return itsPrototype1;
    }

  else if (CAT2 == cat)
    {
      return itsPrototype2;
    }

  else
    throw Util::Error("unknown category enumerator in findStoredExemplar");

  return Mtx::emptyMtx(); // can't happen, but placate the compiler
}

static const char vcid_cmodelwpsm_cc[] = "$Header$";
#endif // !CMODELWPSM_CC_DEFINED
