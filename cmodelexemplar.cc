///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 14:32:31 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_CC_DEFINED
#define CMODELEXEMPLAR_CC_DEFINED

#include "cmodel/cmodelexemplar.h"
#include "cmodel/minkbinder.h"

#include "mtx/mtx.h"
#include "mtx/num.h"

#include "mx/mxwrapper.h"

#include "util/error.h"
#include "util/strings.h"

#include "util/trace.h"

#include <cmath>
#include <vector>

///////////////////////////////////////////////////////////////////////
//
// CModelExemplar member definitions
//
///////////////////////////////////////////////////////////////////////

CModelExemplar::CModelExemplar(const mtx& objParams,
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
                              const mtx& allModelParams,
                              const mx_wrapper& extraArgs)
{
DOTRACE("CModelExemplar::handleRequest");

  if ( action == "getStoredExemplars" )
    {
      mtx category_ = extraArgs.get_field("category").as_mtx();

      int category = int(category_.at(0));

      slice otherParams =
        allModelParams
        .column(0)
        (range(DIM_OBJ_PARAMS+2, allModelParams.ncols()));

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


void CModelExemplar::computeDiffEv(const mtx& objects,
                                   slice& modelParams, mtx& diffEvOut)
{
DOTRACE("CModelExemplar::computeDiffEv");

  const bool newObjects = (objects != itsObjectsCache);

  itsObjectsCache = objects;

  if (newObjects)
    {
      itsEvidence1Cache = mtx(itsNumStoredExemplars, itsObjectsCache.mrows());
      itsEvidence2Cache = mtx(itsNumStoredExemplars, itsObjectsCache.mrows());
    }

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  slice attWeights = modelParams(range(0, DIM_OBJ_PARAMS));

  attWeights.apply(std::abs);

  const bool newAttWts = (itsAttWtsCache.column(0) != attWeights);

  itsAttWtsCache.column(0) = attWeights;

  slice otherParams =
    modelParams(range(DIM_OBJ_PARAMS+2, modelParams.nelems()));

  loadModelParams(otherParams);

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  const mtx& stored1 = getStoredExemplars(CAT1);
  const mtx& stored2 = getStoredExemplars(CAT2);

  mtx_const_iter attWts = attWeights.begin();

  std::vector<mtx_const_iter> exemplars;

  for (int yy = 0; yy < objects.mrows(); ++yy) {
    exemplars.push_back(objects.row_iter(yy));
  }

  for (int x = 0; x < itsNumStoredExemplars; ++x) {
  DOTRACE("minkowski loop");

    bool compute1 = newObjects || newAttWts ||
      (itsStored1Cache.row(x) != stored1.row(x));

    bool compute2 = newObjects || newAttWts ||
      (itsStored2Cache.row(x) != stored2.row(x));

    if (compute1) {
    DOTRACE("compute1");

      const mtx_iter distrust1 = itsEvidence1Cache.row_iter(x);

      MinkowskiBinder binder1(attWts, stored1.row_iter(x),
                              minkPower, minkPowerInv);

      int y = 0;
      for (mtx_iter iter1 = distrust1; iter1.has_more(); ++y, ++iter1) {

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

      const mtx_iter distrust2 = itsEvidence2Cache.row_iter(x);

      MinkowskiBinder binder2(attWts, stored2.row_iter(x),
                              minkPower, minkPowerInv);

      int y = 0;
      for (mtx_iter iter2 = distrust2; iter2.has_more(); ++y, ++iter2) {
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

  diffEvOut.clear(0.0);

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

    const mtx_iter distrust1 = itsEvidence1Cache.row_iter(x);
    const mtx_iter distrust2 = itsEvidence2Cache.row_iter(x);

    const mtx_iter diffEv = diffEvOut.column_iter(0);

    for (mtx_iter iter1 = distrust1, iter2 = distrust2, diff = diffEv;
         iter1.has_more() && diff.has_more();
         ++iter1, ++iter2, ++diff) {
      *diff += *iter1 - *iter2;
    }
  }
}

double CModelExemplar::computeSigmaNoise(double rawSigma) const
{
  return rawSigma * sqrt(itsNumStoredExemplars*2.0);
}

void CModelExemplar::loadModelParams(slice& modelParams) {}

static const char vcid_cmodelexemplar_cc[] = "$Header$";
#endif // !CMODELEXEMPLAR_CC_DEFINED
