///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.h
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 14:31:31 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_H_DEFINED
#define CMODELEXEMPLAR_H_DEFINED

#include "cmodel/classifier.h"

#include "mtx/mtx.h"

class CModelExemplar : public Classifier
{
public:

  enum Category { CAT1, CAT2 };

  enum TransferFunction { EXP_DECAY, LINEAR_DECAY };

  static const int MAX_STORED = -1;

  CModelExemplar(const mtx& objParams,
                 int numStoredExemplars,
                 TransferFunction transferFunc);

  virtual ~CModelExemplar();

  // Handles the request via chain-of-responsibility. Subclasses must
  // be sure to call the superclass version before attempting to
  // process the request.
  virtual RequestResult handleRequest(rutz::fstring action,
                                      const mtx& modelParams,
                                      const mx_wrapper& extraArgs);

  int numStoredExemplars() const { return itsNumStoredExemplars; }

  TransferFunction transferFunction() const { return itsTransferFunc; }

protected:
  const mtx& training1() const { return itsTraining1; }
  const mtx& training2() const { return itsTraining2; }

  int numTrainingExemplars() const { return itsNumTrainingExemplars; }

private:
  const mtx itsTraining1;
  const mtx itsTraining2;

  const int itsNumTrainingExemplars;

  const int itsNumStoredExemplars;

  const TransferFunction itsTransferFunc;

  mtx itsObjectsCache;

  mtx itsStored1Cache;
  mtx itsStored2Cache;

  mtx itsEvidence1Cache;
  mtx itsEvidence2Cache;

  mtx itsAttWtsCache;

  virtual void computeDiffEv(const mtx& objects,
                             slice& modelParams, mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

  virtual void loadModelParams(slice& modelParams);

  virtual const mtx& getStoredExemplars(Category cat) = 0;
};

static const char vcid_cmodelexemplar_h[] = "$Header$";
#endif // !CMODELEXEMPLAR_H_DEFINED
