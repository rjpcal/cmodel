///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Mon Mar 12 16:29:27 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

class Mtx;

class ConstSlice;

class Classifier {
private:
  const Mtx& itsObjParams;
  int itsNumAllExemplars;
  const Mtx& itsObservedIncidence;
  double* const itsDiffEvidence;
  double* const itsPredictedProbability;

  void forwardProbit(double thresh, double sigmaNoise) const;

  enum LogLType { CURRENT, FULL };

  double computeLogL(LogLType type);

  void resetDiffEv()
  {
    for (int i = 0; i < itsNumAllExemplars; ++i)
      itsDiffEvidence[i] = 0.0;
  }

  // Must be overridden by subclasses
  virtual void computeDiffEv(Mtx& modelParams) = 0;

  virtual double fetchSigmaNoise(const Mtx& modelParams) const = 0;

public:
  Classifier(const Mtx& objParams, const Mtx& observedIncidence);
  virtual ~Classifier();

  double currentLogL(Mtx& modelParams);  

  double fullLogL();

  double deviance(Mtx& modelParams);

protected:
  double& diffEvidence(int i) { return itsDiffEvidence[i]; }
  int numAllExemplars() const { return itsNumAllExemplars; }

  ConstSlice exemplar(int i) const;
  int exemplarCategory(int i) const;

  static const int DIM_OBJ_PARAMS = 4;
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
