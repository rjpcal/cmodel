///////////////////////////////////////////////////////////////////////
//
// cmodelrxm.cc
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Jul 31 14:47:40 2002
// written: Wed Jul 31 14:47:40 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELRXM_CC_DEFINED
#define CMODELRXM_CC_DEFINED

#include "cmodel/cmodelrxm.h"

#include "cmodel/cmodelutil.h"

#include "util/error.h"

#include "mtx/mtx.h"

#include "util/trace.h"

CModelRxm::CModelRxm(const Mtx& objParams,
		     TransferFunction transferFunc,
		     int numStoredExemplars) :
  CModelExemplar(objParams, numStoredExemplars, transferFunc),
  itsStored1(numStoredExemplars, DIM_OBJ_PARAMS),
  itsStored2(numStoredExemplars, DIM_OBJ_PARAMS),
  itsHiLo1(CModelUtil::getHiLo(objectsOfCategory(0))),
  itsHiLo2(CModelUtil::getHiLo(objectsOfCategory(1)))
{
DOTRACE("CModelRxm::CModelRxm");
}

CModelRxm::~CModelRxm() {}

void CModelRxm::loadModelParams(Slice& modelParams)
{
DOTRACE("CModelRxm::loadModelParams");

  const int nex = numStoredExemplars();

  const Mtx storedExemplars =
    CModelUtil::getStoredExemplars(modelParams, nex, itsHiLo1, itsHiLo2);

  itsStored1 = storedExemplars.sub(row_range_n(0, nex));
  itsStored2 = storedExemplars.sub(row_range_n(nex, nex));
}

const Mtx& CModelRxm::getStoredExemplars(Category cat)
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
    throw Util::Error("unknown category enumerator in findStoredExemplar");

  return Mtx::emptyMtx(); // can't happen, but placate the compiler
}

static const char vcid_cmodelrxm_cc[] = "$Header$";
#endif // !CMODELRXM_CC_DEFINED
