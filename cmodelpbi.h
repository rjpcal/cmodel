///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 18:08:46 2001
// written: Wed Feb 20 18:09:40 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_H_DEFINED
#define CMODELPBI_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(CLASSIFIER_H_DEFINED)
#include "classifier.h"
#endif

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
