///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.cc
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:09:09 2001
// written: Mon Feb  4 14:01:24 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_CC_DEFINED
#define CMODELPBI_CC_DEFINED

#include "cmodelpbi.h"

#include "mtx.h"

#include "util/trace.h"

CModelPbi::CModelPbi(const Mtx& objParams) :
  Classifier(objParams)
{}

CModelPbi::~CModelPbi() {}

void CModelPbi::computeDiffEv(const Mtx& objects,
                              Slice& modelParams, Mtx& diffEvOut)
{
DOTRACE("CModelPbi::computeDiffEv");

  MtxConstIter attWeights =
    static_cast<const Slice&>(modelParams.leftmost(DIM_OBJ_PARAMS)).begin();

  int i = 0;
  MtxIter diffEv = diffEvOut.columnIter(0);

  for (; diffEv.hasMore(); ++i, ++diffEv)
    *diffEv = -1.0 * innerProduct(attWeights, objects.rowIter(i));
}

double CModelPbi::computeSigmaNoise(double /* rawSigma */) const
{
  return 1.0;
}

static const char vcid_cmodelpbi_cc[] = "$Header$";
#endif // !CMODELPBI_CC_DEFINED
