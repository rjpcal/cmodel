///////////////////////////////////////////////////////////////////////
//
// cmodelcuevalidity.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Tue Apr 10 09:47:00 2001
// written: Wed Jul 31 15:06:39 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCUEVALIDITY_H_DEFINED
#define CMODELCUEVALIDITY_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(CLASSIFIER_H_DEFINED)
#include "cmodel/classifier.h"
#endif

class CModelCueValidity : public Classifier
{
public:
  enum Flag { FREQ_WEIGHT, NO_FREQ_WEIGHT };

  CModelCueValidity(const Mtx& objParams, Flag f);

  virtual ~CModelCueValidity();

protected:
  virtual void computeDiffEv(const Mtx& objects,
                             Slice& modelParams, Mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

private:
  Flag itsFlags;
  const Mtx itsTraining1;
  const Mtx itsTraining2;
};

static const char vcid_cmodelcuevalidity_h[] = "$Header$";
#endif // !CMODELCUEVALIDITY_H_DEFINED
