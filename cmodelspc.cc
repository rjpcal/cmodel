///////////////////////////////////////////////////////////////////////
//
// cmodelspc.cc
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Mon Feb  4 14:01:03 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELSPC_CC_DEFINED
#define CMODELSPC_CC_DEFINED

#include "cmodel/cmodelspc.h"

#include "cmodel/cmodelutil.h"
#include "cmodel/eucbinder.h"

#include "util/error.h"
#include "util/strings.h"

#include <cmath>
#include <limits>

#include "util/trace.h"

CModelSPC::CModelSPC(const Mtx& objParams, int numStoredExemplars) :
  Classifier(objParams),
  itsNumStoredExemplars(numStoredExemplars),
  itsHiLo0(CModelUtil::getHiLo(objectsOfCategory(0))),
  itsHiLo1(CModelUtil::getHiLo(objectsOfCategory(1)))
{}

CModelSPC::~CModelSPC() {}

int CModelSPC::numModelParams() const
{
DOTRACE("CModelSPC::numModelParams");

  return Classifier::numModelParams()
    + (itsNumStoredExemplars * 2 * DIM_OBJ_PARAMS);
}

Classifier::RequestResult
CModelSPC::handleRequest(fstring action,
                         const Mtx& allModelParams,
                         const MxWrapper& extraArgs)
{
DOTRACE("CmodelSPC::handleRequest");

  if ( action == "getStoredExemplars" )
    {
      Mtx category_ = extraArgs.getField("category").asMtx();

      int category = int(category_.at(0));

      slice modelParams = allModelParams.column(0);

      Mtx storedExemplars =
        CModelUtil::getStoredExemplars(modelParams,
                                       itsNumStoredExemplars,
                                       itsHiLo0,
                                       itsHiLo1);

      if (category == 0)
        return storedExemplars(row_range_n(0, itsNumStoredExemplars));

      if (category == 1)
        return storedExemplars(row_range_n(itsNumStoredExemplars,
                                           itsNumStoredExemplars));

      throw Util::Error("unknown category while processing request"
                        "'getStoredExemplars'");
    }

  return Classifier::handleRequest(action, allModelParams, extraArgs);
}

void CModelSPC::computeDiffEv(const Mtx& objects,
                              slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelSPC::computeDiffEv");

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  Mtx attWeights(modelParams(range(0, DIM_OBJ_PARAMS)));

  attWeights.apply(std::abs);

  const slice otherParams =
    modelParams(range(Classifier::DIM_OBJ_PARAMS+2,
                      modelParams.nelems()));

  const Mtx storedExemplars =
    CModelUtil::getStoredExemplars(otherParams,
                                   itsNumStoredExemplars,
                                   itsHiLo0,
                                   itsHiLo1);

  {DOTRACE("<cmodelspc.cc>::loop");

  // Loop over the test objects
  for (int r = 0; r < objects.mrows(); ++r)
    {
      EuclideanBinder<DIM_OBJ_PARAMS>
        ebinder(mtx_const_iter(attWeights.columnIter(0)), objects.rowIter(r));

      // compute distances between test object and stored exemplars, looking
      // for the nearest exemplar from each category

      double mindist0 = std::numeric_limits<double>::max();
      double mindist1 = std::numeric_limits<double>::max();

      {DOTRACE("<cmodelspc.cc>::minkDist2");

      for (int rr = 0, rr2 = itsNumStoredExemplars;
           rr < itsNumStoredExemplars;
           ++rr, ++rr2)
        {
          double dist0 = ebinder.eucDist(storedExemplars.rowIter(rr));
          double dist1 = ebinder.eucDist(storedExemplars.rowIter(rr2));

          if (dist0 < mindist0) mindist0 = dist0;
          if (dist1 < mindist1) mindist1 = dist1;
        }
      }

      // the evidence is the difference of the distances to the nearest
      // exemplars of each category
      diffEvOut.at(r) = (mindist1 - mindist0);
    }
  }
}

double CModelSPC::computeSigmaNoise(double rawSigma) const
{
  return rawSigma;
}

int CModelSPC::fillModelParamsBounds(Mtx& bounds, int startRow) const
{
DOTRACE("CModelSPC::fillModelParamsBounds");

  startRow += Classifier::fillModelParamsBounds(bounds, startRow);

  // number per category
  const int n = itsNumStoredExemplars * DIM_OBJ_PARAMS;

  const int start1 = startRow;
  const int start2 = startRow+n;

  for (int row = 0; row < n; row += DIM_OBJ_PARAMS)
    {
      // first category lower bound
      bounds.sub(row_range_n(start1+row, DIM_OBJ_PARAMS),
                 col_range_n(0, 1))
        =
        itsHiLo0.sub(row_range_n(0,1));

      // first category upper bound
      bounds.sub(row_range_n(start1+row, DIM_OBJ_PARAMS),
                 col_range_n(1, 1))
        =
        itsHiLo0.sub(row_range_n(1,1));

      // second category lower bound
      bounds.sub(row_range_n(start2+row, DIM_OBJ_PARAMS),
                 col_range_n(0, 1))
        =
        itsHiLo1.sub(row_range_n(0,1));

      // second category upper bound
      bounds.sub(row_range_n(start2+row, DIM_OBJ_PARAMS),
                 col_range_n(1, 1))
        =
        itsHiLo1.sub(row_range_n(1,1));
    }

  return startRow + n*2;
}

static const char vcid_cmodelspc_cc[] = "$Header$";
#endif // !CMODELSPC_CC_DEFINED
