///////////////////////////////////////////////////////////////////////
//
// cmodelcuevalidity.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Tue Apr 10 09:47:56 2001
// written: Tue Apr 10 11:20:07 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCUEVALIDITY_CC_DEFINED
#define CMODELCUEVALIDITY_CC_DEFINED

#include "cmodelcuevalidity.h"

#include "mtx.h"

#include "trace.h"

#include "mex.h"

namespace {

  int countCue(double cue, MtxConstIter itr)
  {
	 int result = 0;
	 while (itr.hasMore())
		{
		  if (*itr == cue) ++result;
		  ++itr;
		}
	 return result;
  }
}

CModelCueValidity::CModelCueValidity(const Mtx& objParams,
												 const Mtx& observedIncidence,
												 Flag f) :
  Classifier(objParams, observedIncidence),
  itsFlags(f),
  itsTraining1(objectsOfCategory(0)),
  itsTraining2(objectsOfCategory(1))
{}

CModelCueValidity::~CModelCueValidity() {}

void CModelCueValidity::computeDiffEv(const Mtx& objects,
												  Slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelCueValidity::computeDiffEv");

  Mtx pCond(2,DIM_OBJ_PARAMS);

  Mtx attWeights(modelParams.leftmost(DIM_OBJ_PARAMS));
  attWeights.apply(abs);

  double nTrainers = itsTraining1.mrows() + itsTraining2.mrows();

  for (int i = 0; i < objects.mrows(); ++i)
	 {
		const Slice currentObj = objects.row(i);

		for (int d = 0; d < DIM_OBJ_PARAMS; ++d)
		  {
			 int count1 = countCue(currentObj[d], itsTraining1.columnIter(d));
			 int count2 = countCue(currentObj[d], itsTraining2.columnIter(d));

			 double pJoint1 = double(count1) / nTrainers;
			 double pJoint2 = double(count2) / nTrainers;

			 double pPrior =
				countCue(currentObj[d], objects.columnIter(d)) / nTrainers;

			 double pCond1 = pJoint1 / pPrior;
			 double pCond2 = pJoint2 / pPrior;

			 double weight = (itsFlags == NO_FREQ_WEIGHT) ?
				0.0 :
				1.0 / (1.0+count1+count2);

			 pCond.at(0,d) = 0.5*weight + pCond1*(1.0-weight);
			 pCond.at(1,d) = 0.5*weight + pCond2*(1.0-weight);
		  }

		Mtx cueValidity(2,1); cueValidity.assign_MMmul(pCond, attWeights);

		diffEvOut.at(i,0) = cueValidity.at(0,0) - cueValidity.at(1,0);
	 }
}

double CModelCueValidity::computeSigmaNoise(double rawSigma) const
{
  return rawSigma * sqrt(2.0);
}

static const char vcid_cmodelcuevalidity_cc[] = "$Header$";
#endif // !CMODELCUEVALIDITY_CC_DEFINED
