///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:08:46 2001
// written: Fri Mar  9 18:29:54 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_H_DEFINED
#define CMODELPBI_H_DEFINED

#include "classifier.h"

class CModelPbi : public Classifier {
public:
  CModelPbi(const Rat& objParams, const Rat& observedIncidence);

  virtual ~CModelPbi();

private:
  virtual void computeDiffEv(Rat& modelParams);
  virtual double fetchSigmaNoise(const Rat& modelParams) const;
};

static const char vcid_cmodelpbi_h[] = "$Header$";
#endif // !CMODELPBI_H_DEFINED
