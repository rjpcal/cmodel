///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:08:46 2001
// written: Fri Apr  6 10:27:20 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_H_DEFINED
#define CMODELPBI_H_DEFINED

#include "classifier.h"

class CModelPbi : public Classifier {
public:
  CModelPbi(const Mtx& objParams, const Mtx& observedIncidence);

  virtual ~CModelPbi();

private:
  virtual void computeDiffEv(Slice& modelParams, Mtx& diffEvOut);
  virtual double computeSigmaNoise(double rawSigma) const;
};

static const char vcid_cmodelpbi_h[] = "$Header$";
#endif // !CMODELPBI_H_DEFINED
