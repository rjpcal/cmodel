///////////////////////////////////////////////////////////////////////
//
// cmodelcuevalidity.cc
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Tue Apr 10 09:47:56 2001
// written: Mon Mar  4 12:02:04 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCUEVALIDITY_CC_DEFINED
#define CMODELCUEVALIDITY_CC_DEFINED

#include "cmodelcuevalidity.h"

#include "mtx/mtx.h"

#include "util/trace.h"

#include <map>

namespace
{
  void countCue(std::map<double, int>& counts, MtxConstIter src)
  {
    while (src.hasMore())
      {
        counts[*src] += 1;
        ++src;
      }
  }
}

CModelCueValidity::CModelCueValidity(const Mtx& objParams,
                                     Flag f) :
  Classifier(objParams),
  itsFlags(f),
  itsTraining1(objectsOfCategory(0)),
  itsTraining2(objectsOfCategory(1))
{}

CModelCueValidity::~CModelCueValidity() {}

void CModelCueValidity::computeDiffEv(const Mtx& objects,
                                      Slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelCueValidity::computeDiffEv");

  Mtx attWeights(modelParams(range(0, DIM_OBJ_PARAMS)));
  attWeights.apply(std::abs);

  double nTrainers = itsTraining1.mrows() + itsTraining2.mrows();

  diffEvOut.setAll(0.0);

  for (int d = 0; d < DIM_OBJ_PARAMS; ++d)
    {
      const MtxConstIter objectsColumn = objects.columnIter(d);

      std::map<double, int> t1counts, t2counts, allcounts;
      countCue(t1counts, itsTraining1.columnIter(d));
      countCue(t2counts, itsTraining2.columnIter(d));
      countCue(allcounts, objects.columnIter(d));

      MtxConstIter objectIter = objectsColumn;

      const double attWeight = attWeights.at(d,0);

      MtxIter diffEvIter = diffEvOut.columnIter(0);

      for (; objectIter.hasMore(); ++objectIter, ++diffEvIter)
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
