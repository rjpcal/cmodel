///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:31:31 2001
// written: Tue Mar 13 12:50:30 2001
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

  const Mtx& training1() const { return itsTraining1; }
  const Mtx& training2() const { return itsTraining2; }

  int numTrainingExemplars() const { return itsNumTrainingExemplars; }

private:
  const int itsNumTrainingExemplars;

  Mtx itsTraining1;
  Mtx itsTraining2;

  const int itsNumStoredExemplars;

  void doDiffEvidence(const ConstSlice& attWeights,
							 const ConstSlice& storedExemplar1,
							 const ConstSlice& storedExemplar2,
							 double minkPower,
							 double minkPowerInv);

  void doDiffEvidence2(const ConstSlice& attWeights,
							  const ConstSlice& storedExemplar1,
							  const ConstSlice& storedExemplar2);

  virtual void computeDiffEv(Slice& modelParams);
  virtual double computeSigmaNoise(double rawSigma) const;

  virtual void loadModelParams(Slice& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual ConstSlice findStoredExemplar(Category cat, int n) = 0;
};

static const char vcid_cmodelexemplar_h[] = "$Header$";
#endif // !CMODELEXEMPLAR_H_DEFINED
