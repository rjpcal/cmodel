///////////////////////////////////////////////////////////////////////
//
// cmodelcuevalidity.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Tue Apr 10 09:47:56 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCUEVALIDITY_CC_DEFINED
#define CMODELCUEVALIDITY_CC_DEFINED

#include "cmodel/cmodelcuevalidity.h"

#include "mtx/mtx.h"

#include "util/trace.h"

#include <map>

namespace
{
  void countCue(std::map<double, int>& counts, mtx_const_iter src)
  {
    while (src.has_more())
      {
        counts[*src] += 1;
        ++src;
      }
  }
}

CModelCueValidity::CModelCueValidity(const mtx& objParams,
                                     Flag f) :
  Classifier(objParams),
  itsFlags(f),
  itsTraining1(objectsOfCategory(0)),
  itsTraining2(objectsOfCategory(1))
{}

CModelCueValidity::~CModelCueValidity() {}

void CModelCueValidity::computeDiffEv(const mtx& objects,
                                      slice& modelParams, mtx& diffEvOut)
{
DOTRACE("CModelCueValidity::computeDiffEv");

  mtx attWeights(modelParams(range(0, DIM_OBJ_PARAMS)));
  attWeights.apply(std::abs);

  double nTrainers = itsTraining1.mrows() + itsTraining2.mrows();

  diffEvOut.clear(0.0);

  for (int d = 0; d < DIM_OBJ_PARAMS; ++d)
    {
      const mtx_const_iter objectsColumn = objects.column_iter(d);

      std::map<double, int> t1counts, t2counts, allcounts;
      countCue(t1counts, itsTraining1.column_iter(d));
      countCue(t2counts, itsTraining2.column_iter(d));
      countCue(allcounts, objects.column_iter(d));

      mtx_const_iter objectIter = objectsColumn;

      const double attWeight = attWeights.at(d,0);

      mtx_iter diffEvIter = diffEvOut.column_iter(0);

      for (; objectIter.has_more(); ++objectIter, ++diffEvIter)
        {
          const int count1 = t1counts[*objectIter];
          const int count2 = t2counts[*objectIter];

          const double pJoint1 = double(count1) / nTrainers;
          const double pJoint2 = double(count2) / nTrainers;

          const double pPrior =
            double(allcounts[*objectIter]) / nTrainers;

          const double weight = (itsFlags == NO_FREQ_WEIGHT) ?
            0.0 :
            1.0 / (1.0+count1+count2);

          const double pCond1 = 0.5*weight + (1.0-weight) * pJoint1 / pPrior;
          const double pCond2 = 0.5*weight + (1.0-weight) * pJoint2 / pPrior;

          *diffEvIter += (attWeight * pCond1) - (attWeight * pCond2);
        }
    }
}

double CModelCueValidity::computeSigmaNoise(double rawSigma) const
{
  return rawSigma * sqrt(2.0);
}

static const char vcid_cmodelcuevalidity_cc[] = "$Header$";
#endif // !CMODELCUEVALIDITY_CC_DEFINED
