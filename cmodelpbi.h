///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 1998-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:08:46 2001
// written: Mon Feb  4 18:12:31 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_H_DEFINED
#define CMODELPBI_H_DEFINED

#include "classifier.h"

class CModelPbi : public Classifier
{
public:
  CModelPbi(const Mtx& objParams);

  virtual ~CModelPbi();

private:
  virtual void computeDiffEv(const Mtx& objects,
                             Slice& modelParams, Mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;
};

static const char vcid_cmodelpbi_h[] = "$Header$";
#endif // !CMODELPBI_H_DEFINED
