///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.cc
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 18:09:09 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_CC_DEFINED
#define CMODELPBI_CC_DEFINED

#include "cmodel/cmodelpbi.h"

#include "mtx/mtx.h"

#include "util/trace.h"

CModelPbi::CModelPbi(const mtx& objParams) :
  Classifier(objParams)
{}

CModelPbi::~CModelPbi() {}

void CModelPbi::computeDiffEv(const mtx& objects,
                              slice& modelParams, mtx& diffEvOut)
{
DOTRACE("CModelPbi::computeDiffEv");

  mtx_const_iter attWeights =
    static_cast<const slice&>(modelParams(range(0, DIM_OBJ_PARAMS))).begin();

  int i = 0;
  mtx_iter diffEv = diffEvOut.column_iter(0);

  for (; diffEv.has_more(); ++i, ++diffEv)
    *diffEv = -1.0 * inner_product(attWeights, objects.row_iter(i));
}

double CModelPbi::computeSigmaNoise(double /* rawSigma */) const
{
  return 1.0;
}

static const char vcid_cmodelpbi_cc[] = "$Header$";
#endif // !CMODELPBI_CC_DEFINED
