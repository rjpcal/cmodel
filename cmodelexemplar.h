///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:31:31 2001
// written: Fri Mar 16 17:48:37 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_H_DEFINED
#define CMODELEXEMPLAR_H_DEFINED

#include "classifier.h"

#include "mtx.h"

class CModelExemplar : public Classifier {
public:

  enum Category { CAT1, CAT2 };

  enum TransferFunction { EXP_DECAY, LINEAR_DECAY };

  CModelExemplar(const Mtx& objParams,
					  const Mtx& observedIncidence,
					  int numStoredExemplars,
					  TransferFunction transferFunc);

  virtual ~CModelExemplar();

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

  const TransferFunction itsTransferFunc;

  virtual void computeDiffEv(Slice& modelParams);
  virtual double computeSigmaNoise(double rawSigma) const;

  virtual void loadModelParams(Slice& modelParams);

  virtual const Mtx& getStoredExemplars(Category cat) = 0;
};

static const char vcid_cmodelexemplar_h[] = "$Header$";
#endif // !CMODELEXEMPLAR_H_DEFINED
