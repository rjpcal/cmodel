///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Mon Feb 25 13:54:31 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(MTX_H_DEFINED)
#include "mtx/mtx.h"
#endif

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(MXWRAPPER_H_DEFINED)
#include "mx/mxwrapper.h"
#endif

class fstring;

class Classifier
{
public:
  Classifier(const Mtx& objParams);
  virtual ~Classifier();

  // Returns the classification probability for each of 'testObjects'
  Mtx classifyObjects(Slice& modelParams, const Mtx& testObjects);

  double currentLogL(Slice& modelParams,
                     const Mtx& testObjects,
                     const Mtx& observedIncidence);

  double fullLogL(const Mtx& observedIncidence);

  double deviance(Slice& modelParams,
                  const Mtx& testObjects,
                  const Mtx& observedIncidence);

  struct RequestResult {
    RequestResult() : requestHandled(false), result(Mtx(0,0)) {}

    RequestResult(Mtx res) :
      requestHandled(true), result(res) {}

    RequestResult(MxWrapper res) :
      requestHandled(true), result(res) {}

    const bool requestHandled;
    MxWrapper result;
  };

  // Handles the request via chain-of-responsibility. Subclasses must
  // be sure to call the superclass version before attempting to
  // process the request.
  virtual RequestResult handleRequest(fstring action,
                                      const Mtx& modelParams,
                                      const MxWrapper& extraArgs);

  static const int DIM_OBJ_PARAMS = 4;

protected:
  int numAllExemplars() const { return itsNumAllExemplars; }

  Mtx objectsOfCategory(int category) const;

  // Must be overridden by subclasses
  virtual void computeDiffEv(const Mtx& objects,
                             Slice& modelParams, Mtx& diffEvOut) = 0;

  virtual double computeSigmaNoise(double rawSigma) const = 0;


private:
  const Mtx itsObjCategories;
  const Mtx itsObjects;
  int itsNumAllExemplars;
  Mtx itsDiffEvidence;

  Mtx itsObservedIncidenceCache;
  double itsCachedLogL_1_2;

  static Mtx forwardProbit(const Mtx& diffEv,
                           double thresh, double sigmaNoise);

  double computeLogL(const Mtx& predictedProbability,
                     const Mtx& observedIncidence);

  // Count the number of objects with the given category
  int countCategory(int category) const;

  Slice exemplar(int i) const;
  int exemplarCategory(int i) const;
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
