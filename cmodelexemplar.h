///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:31:31 2001
// written: Fri Apr  6 16:44:00 2001
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

  static const int MAX_STORED = -1;

  CModelExemplar(const Mtx& objParams,
					  const Mtx& observedIncidence,
					  int numStoredExemplars,
					  TransferFunction transferFunc);

  virtual ~CModelExemplar();

protected:
  int numStoredExemplars() const { return itsNumStoredExemplars; }

  const Mtx& training1() const { return itsTraining1; }
  const Mtx& training2() const { return itsTraining2; }

  int numTrainingExemplars() const { return itsNumTrainingExemplars; }

private:
  const Mtx itsTraining1;
  const Mtx itsTraining2;

  const int itsNumTrainingExemplars;

  const int itsNumStoredExemplars;

  const TransferFunction itsTransferFunc;

  Mtx itsObjectsCache;

  Mtx itsStored1Cache;
  Mtx itsStored2Cache;

  Mtx itsEvidence1Cache;
  Mtx itsEvidence2Cache;

  Mtx itsAttWtsCache;

  virtual void computeDiffEv(const Mtx& objects,
									  Slice& modelParams, Mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

  virtual void loadModelParams(Slice& modelParams);

  virtual const Mtx& getStoredExemplars(Category cat) = 0;
};

static const char vcid_cmodelexemplar_h[] = "$Header$";
#endif // !CMODELEXEMPLAR_H_DEFINED
