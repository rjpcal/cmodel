///////////////////////////////////////////////////////////////////////
//
// cmodelcuevalidity.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Tue Apr 10 09:47:00 2001
// written: Wed Feb 20 18:07:56 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCUEVALIDITY_H_DEFINED
#define CMODELCUEVALIDITY_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(CLASSIFIER_H_DEFINED)
#include "classifier.h"
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
  Mtx itsTraining1;
  Mtx itsTraining2;
};

static const char vcid_cmodelcuevalidity_h[] = "$Header$";
#endif // !CMODELCUEVALIDITY_H_DEFINED
