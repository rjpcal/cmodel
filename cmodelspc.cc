///////////////////////////////////////////////////////////////////////
//
// cmodelspc.cc
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb  4 14:01:03 2002
// written: Thu Feb 14 11:58:26 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELSPC_CC_DEFINED
#define CMODELSPC_CC_DEFINED

#include "cmodelspc.h"

#include "eucbinder.h"

#include "util/error.h"
#include "util/strings.h"

#include <cmath>
#include <limits>

#include "util/trace.h"

namespace
{
  void clampRows(Mtx& src, int firstrow, int nrows,
                 const Mtx& hilo)
  {
    DOTRACE("<cmodelspc.cc>::clampRows");

    for (int c = 0; c < src.ncols(); ++c)
      {
        double clo = hilo.at(0, c);
        double chi = hilo.at(1, c);
        for (int r = 0; r < nrows; ++r)
          {
            double v = src.at(firstrow+r,c);
            if (v < clo) src.at(firstrow+r,c) = clo;
            if (v > chi) src.at(firstrow+r,c) = chi;
          }
      }
  }

  Mtx getStoredExemplars(Slice& allModelParams, int nstored,
                         const Mtx& hilo0, const Mtx& hilo1)
  {
    DOTRACE("<cmodelspc.cc>::getStoredExemplars");

    Slice otherParams =
      allModelParams.rightmost(allModelParams.nelems()-
                               (Classifier::DIM_OBJ_PARAMS+2));

    // reshape params into a stored exemplar matrix

    Mtx storedExemplars = Mtx(otherParams);

    storedExemplars.reshape(2*nstored, Classifier::DIM_OBJ_PARAMS);

    clampRows(storedExemplars, 0, nstored, hilo0);

    clampRows(storedExemplars, nstored, nstored, hilo1);

    return storedExemplars;
  }

  Mtx getHiLo(const Mtx& src)
  {
    DOTRACE("<cmodelspc.cc>::getHiLo");
    Mtx result(2, src.ncols());

    for (int i = 0; i < src.ncols(); ++i)
      {
        result.at(0, i) = src.column(i).min();
        result.at(1, i) = src.column(i).max();
      }

    return result;
  }
}

CModelSPC::CModelSPC(const Mtx& objParams, int numStoredExemplars) :
  Classifier(objParams),
  itsNumStoredExemplars(numStoredExemplars),
  itsHiLo0(getHiLo(objectsOfCategory(0))),
  itsHiLo1(getHiLo(objectsOfCategory(1)))
{}

CModelSPC::~CModelSPC() {}

Classifier::RequestResult
CModelSPC::handleRequest(fstring action,
                         const Mtx& allModelParams,
                         const MxWrapper& extraArgs)
{
DOTRACE("CmodelSPC::handleRequest");

  if ( action == "getStoredExemplars" )
    {
      Mtx category_ = extraArgs.getStructField("category").getMtx();

      int category = int(category_.at(0));

      Slice modelParams = allModelParams.column(0);

      Mtx storedExemplars = getStoredExemplars(modelParams,
                                               itsNumStoredExemplars,
                                               itsHiLo0,
                                               itsHiLo1);

      if (category == 0)
        return storedExemplars.rows(0,
                                    itsNumStoredExemplars);

      if (category == 1)
        return storedExemplars.rows(itsNumStoredExemplars,
                                    itsNumStoredExemplars);

      throw Util::Error("unknown category while processing request"
                        "'getStoredExemplars'");
    }

  return Classifier::handleRequest(action, allModelParams, extraArgs);
}

void CModelSPC::computeDiffEv(const Mtx& objects,
                              Slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelSPC::computeDiffEv");

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  Slice attWeights = modelParams.leftmost(DIM_OBJ_PARAMS);

  attWeights.apply(abs);

  Mtx storedExemplars = getStoredExemplars(modelParams,
                                           itsNumStoredExemplars,
                                           itsHiLo0,
                                           itsHiLo1);

  {DOTRACE("<cmodelspc.cc>::loop");

  // Loop over the test objects
  for (int r = 0; r < objects.mrows(); ++r)
    {
      EuclideanBinder<DIM_OBJ_PARAMS>
        ebinder(attWeights.begin(), objects.rowIter(r));

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

static const char vcid_cmodelspc_cc[] = "$Header$";
#endif // !CMODELSPC_CC_DEFINED
