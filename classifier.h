///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Thu Mar  8 16:49:09 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

class Rat;

class fixed_string;
template <class T> class shared_ptr;

class Classifier {
private:
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
  virtual double sigmaScalingFactor() const = 0;

protected:
  double& diffEvidence(int i) { return itsDiffEvidence[i]; }
  int numAllExemplars() const { return itsNumAllExemplars; }

  static const int DIM_OBJ_PARAMS = 4;

public:
  Classifier(const Rat& objParams, const Rat& observedIncidence);
  virtual ~Classifier();

  static shared_ptr<Classifier> make(const fixed_string& whichType,
												 const Rat& objParams,
												 const Rat& observedIncidence,
												 int numStoredExemplars);

  double currentLogL(Rat& modelParams);  

  double fullLogL();

  double deviance(Rat& modelParams);
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
