///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Wed Mar 28 10:24:51 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

#include "mtx.h"

class Classifier {
private:
  const Mtx itsObjParams;
  int itsNumAllExemplars;
  const Mtx itsObservedIncidence;
  Mtx itsDiffEvidence;
  double* const itsPredictedProbability;
  double itsCachedLogL_1_2;

  void forwardProbit(double thresh, double sigmaNoise) const;

  enum LogLType { CURRENT, FULL };

  double computeLogL(LogLType type);

  // Must be overridden by subclasses
  virtual void computeDiffEv(Slice& modelParams, Mtx& diffEvOut) = 0;

  virtual double computeSigmaNoise(double rawSigma) const = 0;

public:
  Classifier(const Mtx& objParams, const Mtx& observedIncidence);
  virtual ~Classifier();

  double currentLogL(Slice& modelParams);  

  double fullLogL();

  double deviance(Slice& modelParams);

protected:
  int numAllExemplars() const { return itsNumAllExemplars; }

  Slice exemplar(int i) const;
  int exemplarCategory(int i) const;

  static const int DIM_OBJ_PARAMS = 4;
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
