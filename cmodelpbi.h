///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:08:46 2001
// written: Tue Mar 13 12:50:30 2001
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
  virtual void computeDiffEv(Slice& modelParams);
  virtual double computeSigmaNoise(double rawSigma) const;
};

static const char vcid_cmodelpbi_h[] = "$Header$";
#endif // !CMODELPBI_H_DEFINED
