///////////////////////////////////////////////////////////////////////
//
// tmex.cc
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Tue Feb 26 09:27:44 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef TMEX_CC_DEFINED
#define TMEX_CC_DEFINED

#include "mx/mexpkg.h"

#include <mex.h>

#include "util/trace.h"

namespace
{
  MexPkg* mexPkg = 0;

  void terminateModule() { delete mexPkg; }

  const int NARGIN = 0;
  const int NARGOUT = 0;
}

void tmex(int nlhs, mxArray* plhs[],
          int nrhs, const mxArray* prhs[])
{
  plhs[0] = mxCreateScalarDouble(mxGetScalar(prhs[0]) + 5.0);
}

extern "C"
void mexFunction(int nlhs, mxArray* plhs[],
                 int nrhs, const mxArray* prhs[])
{
DOTRACE("mexFunction");
  if (mexPkg == 0)
    {
      mexPkg = new MexPkg("tmex", &terminateModule);
      mexPkg->addFcn("tmex", &tmex, NARGIN, NARGOUT);
    }

  return mexPkg->invokeFcn(nlhs, plhs, nrhs, prhs);
}

static const char vcid_tmex_cc[] = "$Id$ $URL$";
#endif // !TMEX_CC_DEFINED
