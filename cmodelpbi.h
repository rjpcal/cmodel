///////////////////////////////////////////////////////////////////////
//
// cmodelpbi.h
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 18:08:46 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELPBI_H_DEFINED
#define CMODELPBI_H_DEFINED

#include "cmodel/classifier.h"

class CModelPbi : public Classifier
{
public:
  CModelPbi(const mtx& objParams);

  virtual ~CModelPbi();

private:
  virtual void computeDiffEv(const mtx& objects,
                             slice& modelParams, mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;
};

static const char vcid_cmodelpbi_h[] = "$Id$ $URL$";
#endif // !CMODELPBI_H_DEFINED
