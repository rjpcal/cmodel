///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:42:01 2001
// written: Tue Oct 30 11:40:24 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_CC_DEFINED
#define CMODELGCM_CC_DEFINED

#include "cmodelgcm.h"

#include "util/error.h"
#include "mtx.h"

#include "util/trace.h"


CModelGcm::CModelGcm(const Mtx& objParams,
                     TransferFunction transferFunc) :
  CModelExemplar(objParams,
                 MAX_STORED,
                 transferFunc)
{
DOTRACE("CModelGcm::CModelGcm");
}


CModelGcm::~CModelGcm() {}


const Mtx& CModelGcm::getStoredExemplars(Category cat)
{
DOTRACE("CModelGcm::findStoredExemplar");

  if (CAT1 == cat)
    {
      return training1();
    }

  else if (CAT2 == cat)
    {
      return training2();
    }

  else
    throw Util::Error("unknown category enumerator in findStoredExemplar");

  return Mtx::emptyMtx(); // can't happen, but placate the compiler
}

static const char vcid_cmodelgcm_cc[] = "$Header$";
#endif // !CMODELGCM_CC_DEFINED
