///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 14:42:01 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_CC_DEFINED
#define CMODELGCM_CC_DEFINED

#include "cmodel/cmodelgcm.h"

#include "mtx/mtx.h"

#include "util/error.h"

#include "util/trace.h"


CModelGcm::CModelGcm(const mtx& objParams,
                     TransferFunction transferFunc) :
  CModelExemplar(objParams,
                 MAX_STORED,
                 transferFunc)
{
DOTRACE("CModelGcm::CModelGcm");
}


CModelGcm::~CModelGcm() {}


const mtx& CModelGcm::getStoredExemplars(Category cat)
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
    throw rutz::error("unknown category enumerator in findStoredExemplar",
                      SRC_POS);

  return mtx::empty_mtx(); // can't happen, but placate the compiler
}

static const char vcid_cmodelgcm_cc[] = "$Header$";
#endif // !CMODELGCM_CC_DEFINED
