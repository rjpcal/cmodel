///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.cc
//
// Copyright (c) 1998-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:32:31 2001
// written: Mon Feb  4 18:12:31 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_CC_DEFINED
#define CMODELEXEMPLAR_CC_DEFINED

#include "cmodelexemplar.h"
#include "minkbinder.h"
#include "mtx.h"
#include "mxwrapper.h"
#include "num.h"

#include "util/error.h"
#include "util/minivec.h"
#include "util/strings.h"

#include "util/trace.h"

#include <cmath>


///////////////////////////////////////////////////////////////////////
//
// CModelExemplar member definitions
//
///////////////////////////////////////////////////////////////////////

CModelExemplar::CModelExemplar(const Mtx& objParams,
                               int numStoredExemplars,
                               TransferFunction transferFunc) :
  Classifier(objParams),
  itsTraining1(objectsOfCategory(0)),
  itsTraining2(objectsOfCategory(1)),
  itsNumTrainingExemplars(itsTraining1.mrows()),
  itsNumStoredExemplars(numStoredExemplars == MAX_STORED ?
                        itsNumTrainingExemplars : numStoredExemplars),
  itsTransferFunc(transferFunc),

  itsObjectsCache(0,0),
  itsStored1Cache(itsNumStoredExemplars, DIM_OBJ_PARAMS),
  itsStored2Cache(itsNumStoredExemplars, DIM_OBJ_PARAMS),
  itsEvidence1Cache(0,0),
  itsEvidence2Cache(0,0),
  itsAttWtsCache(DIM_OBJ_PARAMS,1)
{
  if (itsNumStoredExemplars <= 0)
    throw Util::Error("must have at least one stored exemplar");

  if (itsTraining1.mrows() != itsTraining2.mrows()) {
    throw Util::Error("the two categories must have the "
                      "same number of training exemplars");
  }
}

CModelExemplar::~CModelExemplar()
{}

Classifier::RequestResult
CModelExemplar::handleRequest(fstring action,
                              const Mtx& allModelParams,
                              const MxWrapper& extraArgs)
{
DOTRACE("CModelExemplar::handleRequest");

  if ( action == "getStoredExemplars" )
    {
      Mtx category_ = extraArgs.getStructField("category").getMtx();

      int category = int(category_.at(0));

      Slice modelParams = allModelParams.column(0);

      Slice otherParams = modelParams.rightmost(modelParams.nelems()-
                                                (DIM_OBJ_PARAMS+2));

      loadModelParams(otherParams);

      if (category == 0)
        return getStoredExemplars(CAT1);

      if (category == 1)
        return getStoredExemplars(CAT2);

      throw Util::Error("unknown category while processing request"
                        "'getStoredExemplars'");
    }

  return Classifier::handleRequest(action, allModelParams, extraArgs);
}


void CModelExemplar::computeDiffEv(const Mtx& objects,
                                   Slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelExemplar::computeDiffEv");

  const bool newObjects = (objects != itsObjectsCache);

  itsObjectsCache = objects;

  if (newObjects)
    {
      itsEvidence1Cache = Mtx(itsNumStoredExemplars, itsObjectsCache.mrows());
      itsEvidence2Cache = Mtx(itsNumStoredExemplars, itsObjectsCache.mrows());
    }

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  Slice attWeights = modelParams.leftmost(DIM_OBJ_PARAMS);

  attWeights.apply(abs);

  const bool newAttWts = (itsAttWtsCache.column(0) != attWeights);

  itsAttWtsCache.column(0) = attWeights;

  Slice otherParams = modelParams.rightmost(modelParams.nelems()-
                                            (DIM_OBJ_PARAMS+2));

  loadModelParams(otherParams);

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  const Mtx& stored1 = getStoredExemplars(CAT1);
  const Mtx& stored2 = getStoredExemplars(CAT2);

  MtxConstIter attWts = attWeights.begin();

  minivec<MtxConstIter> exemplars;

  for (int yy = 0; yy < objects.mrows(); ++yy) {
    exemplars.push_back(objects.rowIter(yy));
  }

  for (int x = 0; x < itsNumStoredExemplars; ++x) {
  DOTRACE("minkowski loop");

    bool compute1 = newObjects || newAttWts ||
      (itsStored1Cache.row(x) != stored1.row(x));

    bool compute2 = newObjects || newAttWts ||
      (itsStored2Cache.row(x) != stored2.row(x));

    if (compute1) {
    DOTRACE("compute1");

      const MtxIter distrust1 = itsEvidence1Cache.rowIter(x);

      MinkowskiBinder binder1(attWts, stored1.rowIter(x),
                              minkPower, minkPowerInv);

      int y = 0;
      for (MtxIter iter1 = distrust1; iter1.hasMore(); ++y, ++iter1) {

        if (minkPower == 2.0) {
          if (EXP_DECAY == itsTransferFunc) {
            *iter1 = Num::fastexp7(-binder1.minkDist2(exemplars.at(y)));
          }
          else if (LINEAR_DECAY == itsTransferFunc) {
            *iter1 = -binder1.minkDist2(exemplars.at(y));
          }
        }
        else {
          if (EXP_DECAY == itsTransferFunc) {
            *iter1 = Num::fastexp7(-binder1.minkDist(exemplars.at(y)));
          }
          else if (LINEAR_DECAY == itsTransferFunc) {
            *iter1 = -binder1.minkDist(exemplars.at(y));
          }
        }

      }
      itsStored1Cache.row(x) = stored1.row(x);
    }

    if (compute2) {
    DOTRACE("compute2");
      itsStored2Cache.row(x) = stored2.row(x);

      const MtxIter distrust2 = itsEvidence2Cache.rowIter(x);

      MinkowskiBinder binder2(attWts, stored2.rowIter(x),
                              minkPower, minkPowerInv);

      int y = 0;
      for (MtxIter iter2 = distrust2; iter2.hasMore(); ++y, ++iter2) {
        if (minkPower == 2.0) {
          if (EXP_DECAY == itsTransferFunc) {
            *iter2 = Num::fastexp7(-binder2.minkDist2(exemplars.at(y)));
          }
          else if (LINEAR_DECAY == itsTransferFunc) {
            *iter2 = -binder2.minkDist2(exemplars.at(y));
          }
        }
        else {
          if (EXP_DECAY == itsTransferFunc) {
            *iter2 = Num::fastexp7(-binder2.minkDist(exemplars.at(y)));
          }
          else if (LINEAR_DECAY == itsTransferFunc) {
            *iter2 = -binder2.minkDist(exemplars.at(y));
          }
        }
      }
    }

  }

  diffEvOut.setAll(0.0);

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

    const MtxIter distrust1 = itsEvidence1Cache.rowIter(x);
    const MtxIter distrust2 = itsEvidence2Cache.rowIter(x);

    const MtxIter diffEv = diffEvOut.columnIter(0);

    for (MtxIter iter1 = distrust1, iter2 = distrust2, diff = diffEv;
         iter1.hasMore() && diff.hasMore();
         ++iter1, ++iter2, ++diff) {
      *diff += *iter1 - *iter2;
    }
  }
}

double CModelExemplar::computeSigmaNoise(double rawSigma) const
{
  return rawSigma * sqrt(itsNumStoredExemplars*2.0);
}

void CModelExemplar::loadModelParams(Slice& modelParams) {}

static const char vcid_cmodelexemplar_cc[] = "$Header$";
#endif // !CMODELEXEMPLAR_CC_DEFINED
