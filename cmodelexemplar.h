///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:31:31 2001
// written: Fri Mar  9 16:50:26 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_H_DEFINED
#define CMODELEXEMPLAR_H_DEFINED

#include "classifier.h"

class Slice;

class CModelExemplar : public Classifier {
public:
  CModelExemplar(const Rat& objParams,
					  const Rat& observedIncidence,
					  int numStoredExemplars);

  virtual ~CModelExemplar();

  enum Category { CAT1, CAT2 };

protected:
  int numStoredExemplars() const { return itsNumStoredExemplars; }

  // Count the category training exemplars
  static int countCategory(const Rat& params, int category);

private:
  const Rat& itsObjParams;
  const int itsNumStoredExemplars;

  void doDiffEvidence(const double* attWeights,
							 const Slice& storedExemplar1,
							 const Slice& storedExemplar2,
							 double minkPower,
							 double minkPowerInv);

  void doDiffEvidence2(const double* attWeights,
							  const Slice& storedExemplar1,
							  const Slice& storedExemplar2);

  virtual void computeDiffEv(Rat& modelParams);
  virtual double sigmaScalingFactor() const;

  virtual void loadModelParams(Rat& modelParams) = 0;

  // The result of this function is only valid until the next call to
  // the function
  virtual Slice findStoredExemplar(Category cat, int n) = 0;
};

static const char vcid_cmodelexemplar_h[] = "$Header$";
#endif // !CMODELEXEMPLAR_H_DEFINED
