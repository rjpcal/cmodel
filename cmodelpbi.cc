///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:09:09 2001
// written: Wed Mar 21 14:04:23 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_CC_DEFINED
#define CMODELPBI_CC_DEFINED

#include "cmodelpbi.h"

#include "mtx.h"

#include "trace.h"

CModelPbi::CModelPbi(const Mtx& objParams, const Mtx& observedIncidence) :
  Classifier(objParams, observedIncidence)
{}

CModelPbi::~CModelPbi() {}

void CModelPbi::computeDiffEv(Slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelPbi::computeDiffEv");

  MtxConstIter attWeights =
	 static_cast<const Slice&>(modelParams.leftmost(DIM_OBJ_PARAMS)).begin();

  int i = 0;
  MtxIter diffEv = diffEvOut.colIter(0);

  for (; diffEv.hasMore(); ++i, ++diffEv)
	 *diffEv = -1.0 * innerProduct(attWeights, exemplar(i).begin());
}

double CModelPbi::computeSigmaNoise(double /* rawSigma */) const
{
  return 1.0;
}

static const char vcid_cmodelpbi_cc[] = "$Header$";
#endif // !CMODELPBI_CC_DEFINED
