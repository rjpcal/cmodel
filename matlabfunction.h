///////////////////////////////////////////////////////////////////////
//
// matlabfunction.h
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 18 11:07:28 2002
// written: Mon Feb 18 11:11:58 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MATLABFUNCTION_H_DEFINED
#define MATLABFUNCTION_H_DEFINED

#include "multivarfunction.h"

#include "mtx.h"

#include <libmatlb.h>

///////////////////////////////////////////////////////////////////////
/**
 *
 * MatlabFunction represents a function that takes one variable argument, plus
 * any number of additional arguments that are bound at the time the object is
 * created.
 *
 **/
///////////////////////////////////////////////////////////////////////

class MatlabFunction : public MultivarFunction
{
private:
  mxArray* itsFunfcn;
  mxArray* itsVarargin_ref;

public:
  MatlabFunction(mxArray* funfcn_mx, int nvararg, mxArray** pvararg) :
    MultivarFunction(),
    itsFunfcn(funfcn_mx),
    itsVarargin_ref(0)
  {
    mlfAssign(&itsVarargin_ref, mclCreateVararginCell(nvararg, pvararg));
  }

  virtual ~MatlabFunction()
  {
    mxDestroyArray(itsVarargin_ref);
  }

protected:
  static mxArray* getref(mxArray* varargin)
    {
      return mlfIndexRef(varargin,
                         "{?}",
                         mlfCreateColonIndex());
    }

  double evaluate_mx(mxArray* x_mx);

  virtual double doEvaluate(const Mtx& x)
  {
    return evaluate_mx(x.makeMxArray());
  }
};

static const char vcid_matlabfunction_h[] = "$Header$";
#endif // !MATLABFUNCTION_H_DEFINED
