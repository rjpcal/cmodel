///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Thu Aug  1 11:02:15 2002
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

  /// Returns the number of model params needed by this model.
  /** Subclasses can override to increase the number needed by base
      classes. */
  virtual int numModelParams() const;

  /// Returns an Nx2 matrix of lower and upper bounds on the model params.
  Mtx modelParamsBounds() const;

  /// Returns the classification probability for each of 'testObjects'
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

  /// Handles the request via chain-of-responsibility.
  /** Subclasses must be sure to call the superclass version before
      attempting to process the request. */
  virtual RequestResult handleRequest(fstring action,
                                      const Mtx& modelParams,
                                      const MxWrapper& extraArgs);

  static const int DIM_OBJ_PARAMS = 4;

  const Mtx& objParams() const { return itsObjParams; }

protected:
  int numAllExemplars() const { return itsNumAllExemplars; }

  Mtx objectsOfCategory(int category) const;

  /// Must be overridden by subclasses
  virtual void computeDiffEv(const Mtx& objects,
                             Slice& modelParams, Mtx& diffEvOut) = 0;

  virtual double computeSigmaNoise(double rawSigma) const = 0;

  /// Each class fills in the bounds for its params.
  /** Any subclass that overrides numModelParams() MUST also override this
      function, at least to return a properly adjusted startRow. Subclasses
      must call the base class version first. The return value should be
      the number of params that were filled, so that the caller can adjust
      "startRow" accordingly. The base class version defined here in
      Classifier fills all rows in with extreme lower and upper bounds, so
      if the subclass doesn't need to place any specific bounds on its
      additional parameters, it can simply leave the bounds matrix alone,
      and just return an adjusted startRow. */
  virtual int fillModelParamsBounds(Mtx& bounds, int startRow) const;

private:
  const Mtx itsObjParams; // col 0 --> categories, cols 1-4 --> objects
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
