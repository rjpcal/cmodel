///////////////////////////////////////////////////////////////////////
//
// matlabfunction.h
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Mon Feb 18 11:07:28 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MATLABFUNCTION_H_DEFINED
#define MATLABFUNCTION_H_DEFINED

#include "cmodel/multivarfunction.h"

#include "util/strings.h"

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
  const fstring         itsFuncName;
  const int             itsNvararg;
  const mxArray** const itsPvararg;
  const bool            itsCanUseMatrix;
  mxArray**             itsPrhs;

  Mtx evaluateParallel(const Mtx& models) const;
  Mtx evaluateSerial(const Mtx& models) const;

public:
  MatlabFunction(const fstring& funcName,
                 int nvararg,
                 const mxArray** pvararg,
                 bool canUseMatrix = false);

  virtual ~MatlabFunction();

protected:
  virtual Mtx doEvaluateEach(const Mtx& models);

  virtual double doEvaluate(const Mtx& model);
};

static const char vcid_matlabfunction_h[] = "$Header$";
#endif // !MATLABFUNCTION_H_DEFINED
