///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:08:46 2001
// written: Mon Mar 12 12:31:47 2001
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
  virtual void computeDiffEv(Mtx& modelParams);
  virtual double fetchSigmaNoise(const Mtx& modelParams) const;
};

static const char vcid_cmodelpbi_h[] = "$Header$";
#endif // !CMODELPBI_H_DEFINED
