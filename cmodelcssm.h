///////////////////////////////////////////////////////////////////////
//
// modelcssm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:24:41 2001
// written: Thu Mar  8 16:24:53 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MODELCSSM_H_DEFINED
#define MODELCSSM_H_DEFINED

#include "classifier.h"

class ModelCssm : public Classifier {
private:
  const Rat& itsObjParams;
  const int itsNumStoredExemplars;
  const int itsNum1;
  const int itsNum2;
  const int itsNumTrainingExemplars;

  typedef const double* constDblPtr;

  constDblPtr* const itsCat1;
  constDblPtr* const itsCat2;

  // Count the category training exemplars
  int countCategory(const Rat& params, int category);

  // Scales the weights in place; weights is an input/output argument
  void scaleWeights(double* weights, int numRawWeights);

  void computeSimilarity(const double* attWeights,
								 const double* storedExemplar1,
								 const double* storedExemplar2,
								 double minkPower,
								 double minkPowerInv);

  void computeSimilarity2(const double* attWeights,
								  const double* storedExemplar1,
								  const double* storedExemplar2);

  virtual void computeDiffEv(Rat& modelParams);
  virtual double sigmaScalingFactor() const;

public:
  ModelCssm(const Rat& objParams,
				const Rat& observedIncidence,
				int numStoredExemplars);

  ~ModelCssm();
};

static const char vcid_modelcssm_h[] = "$Header$";
#endif // !MODELCSSM_H_DEFINED
