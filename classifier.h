///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Thu Mar  8 11:54:49 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CLASSIFIER_H_DEFINED
#define CLASSIFIER_H_DEFINED

#ifndef MLF_V2
#define MLF_V2 1
#endif

class Rat;

class Classifier {
public:
  static void forwardProbit(double* diffEvidence,
									 int num_points,
									 double thresh_,
									 double sigmaNoise_,
									 double* ppOut);

  static double loglikelihood(const double* predictedProbability1,
										int numPoints,
										const double* observedIncidence1,
										const double* observedIncidence2);
};

class ModelCssm : public Classifier {
private:
  const Rat& itsObjParams;
  const Rat& itsObservedIncidence;
  const int itsNumStoredExemplars;
  const int itsNum1;
  const int itsNum2;
  const int itsNumTrainingExemplars;
  const int itsNumAllExemplars;

  const int itsDimObjParams;

  double* const itsDiffEvidence;

  typedef const double* constDblPtr;

  constDblPtr* const itsCat1;
  constDblPtr* const itsCat2;

  double* itsPredictedProbability;

  // Count the category training exemplars
  int countCategory(const Rat& params, int category);

  void reset()
  {
    for (int i = 0; i < itsNumAllExemplars; ++i)
      itsDiffEvidence[i] = 0.0;
  }

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

public:
  ModelCssm(const Rat& objParams,
				const Rat& observedIncidence,
				int numStoredExemplars);

  ~ModelCssm();

  double loglikelihoodFor(Rat& modelParams);
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
