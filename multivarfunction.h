///////////////////////////////////////////////////////////////////////
//
// multivarfunction.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 06:45:02 2001
// written: Wed Apr 18 06:48:46 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MULTIVARFUNCTION_H_DEFINED
#define MULTIVARFUNCTION_H_DEFINED

class Mtx;

class MultivarFunction {
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
