///////////////////////////////////////////////////////////////////////
//
// multivarfunction.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Apr 18 06:45:02 2001
// written: Tue Feb 19 09:18:31 2002
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

public:
  MultivarFunction();

  virtual ~MultivarFunction();

  int evalCount() const { return itsEvalCount; }

  double evaluate(const Mtx& x)
  {
    ++itsEvalCount;
    return doEvaluate(x);
  }

  /** Treat each column of x as a set of parameters to be evaluated, and
      return a vector of the results of evaluating each column. This may be
      more efficient than making a number of repeated calls to evaluate(). */
  Mtx evaluateEach(const Mtx& x);

protected:
  virtual double doEvaluate(const Mtx& x) = 0;

  /** May be overridden to provide a more optimized parallel implementation,
      but the default version just makes repeated calls to doEvaluate(). */
  virtual Mtx doEvaluateEach(const Mtx& x);
};

static const char vcid_multivarfunction_h[] = "$Header$";
#endif // !MULTIVARFUNCTION_H_DEFINED
