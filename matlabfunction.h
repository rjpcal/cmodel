///////////////////////////////////////////////////////////////////////
//
// matlabfunction.h
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb 18 11:07:28 2002
// written: Wed Jul 31 15:09:03 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MATLABFUNCTION_H_DEFINED
#define MATLABFUNCTION_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(MULTIVARFUNCTION_H_DEFINED)
#include "cmodel/multivarfunction.h"
#endif

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(STRINGS_H_DEFINED)
#include "util/strings.h"
#endif

typedef struct mxArray_tag mxArray;

///////////////////////////////////////////////////////////////////////
/**
 *
 * MatlabFunction represents a function that takes one variable argument, plus
 * any number of trailing arguments that are bound at the time the object is
 * created.
 *
 **/
///////////////////////////////////////////////////////////////////////

class MatlabFunction : public MultivarFunction
{
private:
  const fstring itsFuncName;
  const int itsNvararg;
  mxArray** const itsPvararg;
  const bool itsCanUseMatrix;

  mxArray** itsPrhs;

  Mtx evaluateParallel(const Mtx& models) const;
  Mtx evaluateSerial(const Mtx& models) const;

public:
  MatlabFunction(const fstring& funcName,
                 int nvararg, mxArray** pvararg,
                 bool canUseMatrix = false);

  virtual ~MatlabFunction();

protected:
  virtual Mtx doEvaluateEach(const Mtx& models);

  virtual double doEvaluate(const Mtx& model);
};

static const char vcid_matlabfunction_h[] = "$Header$";
#endif // !MATLABFUNCTION_H_DEFINED
