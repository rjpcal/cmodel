///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 2001-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Thu Mar  8 09:48:36 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#include "mtx/mtx.h"

#include "mx/mxwrapper.h"

namespace rutz
{
  class fstring;
}

class Classifier
{
public:
  Classifier(const mtx& objParams);
  virtual ~Classifier();

  /// Returns the number of model params needed by this model.
  /** Subclasses can override to increase the number needed by base
      classes. */
  virtual int numModelParams() const;

  /// Returns an Nx2 matrix of lower and upper bounds on the model params.
  mtx modelParamsBounds() const;

  /// Returns the classification probability for each of 'testObjects'
  mtx classifyObjects(slice& modelParams, const mtx& testObjects);

  double currentLogL(slice& modelParams,
                     const mtx& testObjects,
                     const mtx& observedIncidence);

  double fullLogL(const mtx& observedIncidence);

  double deviance(slice& modelParams,
                  const mtx& testObjects,
                  const mtx& observedIncidence);

  struct RequestResult
  {
    RequestResult() :
      requestHandled(false), result(mtx::empty_mtx()) {}

    RequestResult(mtx res) :
      requestHandled(true), result(res) {}

    RequestResult(mx_wrapper res) :
      requestHandled(true), result(res) {}

    const bool requestHandled;
    mx_wrapper result;
  };

  /// Handles the request via chain-of-responsibility.
  /** Subclasses must be sure to call the superclass version before
      attempting to process the request. */
  virtual RequestResult handleRequest(rutz::fstring action,
                                      const mtx& modelParams,
                                      const mx_wrapper& extraArgs);

  static const int DIM_OBJ_PARAMS = 4;

  const mtx& objParams() const { return itsObjParams; }

protected:
  int numAllExemplars() const { return itsNumAllExemplars; }

  mtx objectsOfCategory(int category) const;

  /// Must be overridden by subclasses
  virtual void computeDiffEv(const mtx& objects,
                             slice& modelParams, mtx& diffEvOut) = 0;

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
  virtual int fillModelParamsBounds(mtx& bounds, int startRow) const;

private:
  const mtx itsObjParams; // col 0 --> categories, cols 1-4 --> objects
  const mtx itsObjCategories;
  const mtx itsObjects;
  int itsNumAllExemplars;
  mtx itsDiffEvidence;

  mtx itsObservedIncidenceCache;
  double itsCachedLogL_1_2;

  static mtx forwardProbit(const mtx& diffEv,
                           double thresh, double sigmaNoise);

  double computeLogL(const mtx& predictedProbability,
                     const mtx& observedIncidence);

  // Count the number of objects with the given category
  int countCategory(int category) const;

  slice exemplar(int i) const;
  int exemplarCategory(int i) const;
};

static const char vcid_classifier_h[] = "$Id$ $URL$";
#endif // !CLASSIFIER_H_DEFINED
