///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:09:09 2001
// written: Fri Mar  9 18:32:37 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_CC_DEFINED
#define CMODELPBI_CC_DEFINED

#include "cmodelpbi.h"

#include "rutil.h"

#include "trace.h"

CModelPbi::CModelPbi(const Rat& objParams, const Rat& observedIncidence) :
  Classifier(objParams, observedIncidence)
{}

CModelPbi::~CModelPbi() {}

void CModelPbi::computeDiffEv(Rat& modelParams)
{
DOTRACE("CModelPbi::computeDiffEv");

  Slice attWeights(modelParams.data(), 1);

  for (int i = 0; i < numAllExemplars(); ++i)
	 diffEvidence(i) = -1.0 * Slice::dot(attWeights, exemplar(i),
													 DIM_OBJ_PARAMS);
}

double CModelPbi::fetchSigmaNoise(const Rat& /*modelParams*/) const
{
  return 1.0;
}

static const char vcid_cmodelpbi_cc[] = "$Header$";
#endif // !CMODELPBI_CC_DEFINED
