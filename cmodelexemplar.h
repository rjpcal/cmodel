///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:31:31 2001
// written: Mon Mar 12 12:34:26 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_H_DEFINED
#define CMODELEXEMPLAR_H_DEFINED

#include "classifier.h"

#include "mtx.h"

#include "util/arrays.h"

class CModelExemplar : public Classifier {
public:
  CModelExemplar(const Mtx& objParams,
					  const Mtx& observedIncidence,
					  int numStoredExemplars);

  virtual ~CModelExemplar();

  enum Category { CAT1, CAT2 };

protected:
  int numStoredExemplars() const { return itsNumStoredExemplars; }

  // Count the category training exemplars
  static int countCategory(const Mtx& params, int category);

  const fixed_block<Slice>& training1() const { return itsTraining1; }
  const fixed_block<Slice>& training2() const { return itsTraining2; }

  int numTrainingExemplars() const { return itsNumTrainingExemplars; }

private:
  const int itsNumTrainingExemplars;

  fixed_block<Slice> itsTraining1;
  fixed_block<Slice> itsTraining2;

  const int itsNumStoredExemplars;

  void doDiffEvidence(const double* attWeights,
							 const Slice& storedExemplar1,
							 const Slice& storedExemplar2,
							 double minkPower,
							 double minkPowerInv);

  void doDiffEvidence2(const double* attWeights,
							  const Slice& storedExemplar1,
							  const Slice& storedExemplar2);

  virtual void computeDiffEv(Mtx& modelParams);
  virtual double fetchSigmaNoise(const Mtx& modelParams) const;

  virtual void loadModelParams(Mtx& modelParams) = 0;

  // The result of this function is only valid until the next call to
  // the function
  virtual Slice findStoredExemplar(Category cat, int n) = 0;
};

static const char vcid_cmodelexemplar_h[] = "$Header$";
#endif // !CMODELEXEMPLAR_H_DEFINED
