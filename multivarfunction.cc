///////////////////////////////////////////////////////////////////////
//
// multivarfunction.cc
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Tue Feb 19 09:15:55 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MULTIVARFUNCTION_CC_DEFINED
#define MULTIVARFUNCTION_CC_DEFINED

#include "cmodel/multivarfunction.h"

#include "mtx/mtx.h"

MultivarFunction::MultivarFunction() : itsEvalCount(0) {}

MultivarFunction::~MultivarFunction() {}

Mtx MultivarFunction::evaluateEach(const Mtx& x)
{
  itsEvalCount += x.ncols();

  return doEvaluateEach(x);
}

Mtx MultivarFunction::doEvaluateEach(const Mtx& x)
{
  Mtx result(x.ncols(), 1);

  for (int i = 0; i < x.ncols(); ++i)
    {
      result.at(i) = doEvaluate(Mtx(x.column(i)));
    }

  return result;
}

static const char vcid_multivarfunction_cc[] = "$Header$";
#endif // !MULTIVARFUNCTION_CC_DEFINED
