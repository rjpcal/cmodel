///////////////////////////////////////////////////////////////////////
//
// classifier.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 09:48:36 2001
// written: Thu Mar  8 10:51:19 2001
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
  const Rat& objParams;
  const Rat& observedIncidence;
  const int numStoredExemplars;
  const int num1;
  const int num2;
  const int numTrainingExemplars;
  const int numAllExemplars;

  const int dimObjParams;

  double* const mu;

  typedef const double* constDblPtr;

  constDblPtr* const cat1;
  constDblPtr* const cat2;

  // Count the category training exemplars
  int countCategory(const Rat& params, int category);

  void resetMu()
  {
    for (int i = 0; i < numAllExemplars; ++i)
      mu[i] = 0.0;
  }

  void computeSimilarity(const double* attWeights,
								 const double* storedExemplar1,
								 const double* storedExemplar2,
								 double minkPower,
								 double minkPowerInv);

  void computeSimilarity2(const double* attWeights,
								  const double* storedExemplar1,
								  const double* storedExemplar2);

public:
  ModelCssm(const Rat& objParams_,
				const Rat& observedIncidence_,
				int numStoredExemplars_);

  ~ModelCssm();

  double loglikelihoodFor(Rat& modelParams);
};

static const char vcid_classifier_h[] = "$Header$";
#endif // !CLASSIFIER_H_DEFINED
