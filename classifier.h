///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Fri Mar  9 18:28:23 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

class Rat;
class Slice;

class Classifier {
private:
  const Rat& itsObjParams;
  int itsNumAllExemplars;
  const Rat& itsObservedIncidence;
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
  virtual void computeDiffEv(Rat& modelParams) = 0;

  virtual double fetchSigmaNoise(const Rat& modelParams) const = 0;

public:
  Classifier(const Rat& objParams, const Rat& observedIncidence);
  virtual ~Classifier();

  double currentLogL(Rat& modelParams);  

  double fullLogL();

  double deviance(Rat& modelParams);

protected:
  double& diffEvidence(int i) { return itsDiffEvidence[i]; }
  int numAllExemplars() const { return itsNumAllExemplars; }

  Slice exemplar(int i) const;
  int exemplarCategory(int i) const;

  static const int DIM_OBJ_PARAMS = 4;
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
