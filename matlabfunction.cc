///////////////////////////////////////////////////////////////////////
//
// matlabfunction.cc
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 18 11:10:06 2002
// written: Mon Feb 18 11:10:31 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MATLABFUNCTION_CC_DEFINED
#define MATLABFUNCTION_CC_DEFINED

#include "matlabfunction.h"

#include "util/trace.h"

double MatlabFunction::evaluate_mx(mxArray* x_mx)
{
DOTRACE("MatlabFunction::evaluate_mx");

  mxArray* mx =  mlfFeval(mclValueVarargout(),
                          itsFunfcn,
                          x_mx,
                          getref(itsVarargin_ref),
                          NULL);
  double result = mxGetScalar(mx);
  mxDestroyArray(mx);
  return result;
}

static const char vcid_matlabfunction_cc[] = "$Header$";
#endif // !MATLABFUNCTION_CC_DEFINED
