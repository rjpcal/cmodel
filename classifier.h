///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Fri Apr  6 12:12:56 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#include "mtx.h"

class Classifier {
public:
  Classifier(const Mtx& objParams, const Mtx& observedIncidence);
  virtual ~Classifier();

  // Returns the classification probability for each of 'objects'
  Mtx classifyObjects(Slice& modelParams, const Mtx& objects);

  double currentLogL(Slice& modelParams);  

  double fullLogL();

  double deviance(Slice& modelParams);

protected:
  int numAllExemplars() const { return itsNumAllExemplars; }

  Mtx objectsOfCategory(int category) const;

  static const int DIM_OBJ_PARAMS = 4;

  // Must be overridden by subclasses
  virtual void computeDiffEv(const Mtx& objects,
									  Slice& modelParams, Mtx& diffEvOut) = 0;

  virtual double computeSigmaNoise(double rawSigma) const = 0;


private:
  const Mtx itsObjCategories;
  const Mtx itsObjects;
  int itsNumAllExemplars;
  const Mtx itsObservedIncidence;
  Mtx itsDiffEvidence;
  double itsCachedLogL_1_2;

  Mtx forwardProbit(double thresh, double sigmaNoise) const;

  double computeLogL(const Mtx& predictedProbability);

  // Count the number of objects with the given category
  int countCategory(int category) const;

  Slice exemplar(int i) const;
  int exemplarCategory(int i) const;
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
