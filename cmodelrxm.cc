///////////////////////////////////////////////////////////////////////
//
// cmodelrxm.cc
//
// Copyright (c) 2002-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Wed Jul 31 14:47:40 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELRXM_CC_DEFINED
#define CMODELRXM_CC_DEFINED

#include "cmodel/cmodelrxm.h"

#include "cmodel/cmodelutil.h"

#include "rutz/error.h"

#include "mtx/mtx.h"

#define LOCAL_DEBUG
#include "rutz/debug.h"
#include "rutz/trace.h"

CModelRxm::CModelRxm(const mtx& objParams,
                     TransferFunction transferFunc,
                     int numStoredExemplars) :
  CModelExemplar(objParams, numStoredExemplars, transferFunc),
  itsStored1(mtx::zeros(numStoredExemplars, DIM_OBJ_PARAMS)),
  itsStored2(mtx::zeros(numStoredExemplars, DIM_OBJ_PARAMS)),
  itsHiLo1(CModelUtil::getHiLo(objectsOfCategory(0))),
  itsHiLo2(CModelUtil::getHiLo(objectsOfCategory(1)))
{
GVX_TRACE("CModelRxm::CModelRxm");
}

CModelRxm::~CModelRxm() {}

int CModelRxm::numModelParams() const
{
GVX_TRACE("CModelRxm::numModelParams");

  return CModelExemplar::numModelParams()
    + (numStoredExemplars() * 2 * DIM_OBJ_PARAMS);
}

void CModelRxm::loadModelParams(slice& modelParams)
{
GVX_TRACE("CModelRxm::loadModelParams");

  const int nex = numStoredExemplars();

  const mtx storedExemplars =
    CModelUtil::getStoredExemplars(modelParams, nex, itsHiLo1, itsHiLo2);

  itsStored1 = storedExemplars.sub(row_range_n(0, nex));
  itsStored2 = storedExemplars.sub(row_range_n(nex, nex));
}

const mtx& CModelRxm::getStoredExemplars(Category cat)
{
GVX_TRACE("CModelRxm::getStoredExemplars");
  if (CAT1 == cat)
    {
      return itsStored1;
    }

  else if (CAT2 == cat)
    {
      return itsStored2;
    }

  else
    throw rutz::error("unknown category enumerator in findStoredExemplar",
                      SRC_POS);

  return mtx::empty_mtx(); // can't happen, but placate the compiler
}

int CModelRxm::fillModelParamsBounds(mtx& bounds, int startRow) const
{
GVX_TRACE("CModelRxm::fillModelParamsBounds");

  startRow += CModelExemplar::fillModelParamsBounds(bounds, startRow);

  // number per category
  const int n = numStoredExemplars() * DIM_OBJ_PARAMS;

  const int start1 = startRow;
  const int start2 = startRow+n;

  for (int row = 0; row < n; row += DIM_OBJ_PARAMS)
    {
      // first category lower bound
      bounds.sub(row_range_n(start1+row, DIM_OBJ_PARAMS),
                 col_range_n(0, 1))
        =
        itsHiLo1.sub(row_range_n(0,1));

      // first category upper bound
      bounds.sub(row_range_n(start1+row, DIM_OBJ_PARAMS),
                 col_range_n(1, 1))
        =
        itsHiLo1.sub(row_range_n(1,1));

      // second category lower bound
      bounds.sub(row_range_n(start2+row, DIM_OBJ_PARAMS),
                 col_range_n(0, 1))
        =
        itsHiLo2.sub(row_range_n(0,1));

      // second category upper bound
      bounds.sub(row_range_n(start2+row, DIM_OBJ_PARAMS),
                 col_range_n(1, 1))
        =
        itsHiLo2.sub(row_range_n(1,1));
    }

  return startRow + n*2;
}

static const char vcid_cmodelrxm_cc[] = "$Id$ $URL$";
#endif // !CMODELRXM_CC_DEFINED
