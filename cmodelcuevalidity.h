///////////////////////////////////////////////////////////////////////
//
// cmodelcuevalidity.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Tue Apr 10 09:47:00 2001
// written: Tue Apr 10 11:18:32 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCUEVALIDITY_H_DEFINED
#define CMODELCUEVALIDITY_H_DEFINED

#include "classifier.h"

class CModelCueValidity : public Classifier {
public:
  enum Flag { FREQ_WEIGHT, NO_FREQ_WEIGHT };

  CModelCueValidity(const Mtx& objParams, const Mtx& observedIncidence,
						  Flag f);

  virtual ~CModelCueValidity();

protected:
  virtual void computeDiffEv(const Mtx& objects,
									  Slice& modelParams, Mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

private:
  Flag itsFlags;
  Mtx itsTraining1;
  Mtx itsTraining2;
};

static const char vcid_cmodelcuevalidity_h[] = "$Header$";
#endif // !CMODELCUEVALIDITY_H_DEFINED
