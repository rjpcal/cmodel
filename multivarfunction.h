///////////////////////////////////////////////////////////////////////
//
// multivarfunction.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 06:45:02 2001
// written: Thu Feb 14 11:56:24 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MULTIVARFUNCTION_H_DEFINED
#define MULTIVARFUNCTION_H_DEFINED

class Mtx;

class MultivarFunction
{
private:
  int itsEvalCount;

protected:
  virtual double doEvaluate(const Mtx& x) = 0;

public:
  MultivarFunction() : itsEvalCount(0) {}

  virtual ~MultivarFunction() {}

  int evalCount() const { return itsEvalCount; }

  double evaluate(const Mtx& x)
  {
         ++itsEvalCount;
         return doEvaluate(x);
  }
};

static const char vcid_multivarfunction_h[] = "$Header$";
#endif // !MULTIVARFUNCTION_H_DEFINED
