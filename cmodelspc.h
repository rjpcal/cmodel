///////////////////////////////////////////////////////////////////////
//
// cmodelspc.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Mon Feb  4 13:59:00 2002
// written: Wed Feb 20 18:09:49 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELSPC_H_DEFINED
#define CMODELSPC_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(CLASSIFIER_H_DEFINED)
#include "classifier.h"
#endif

class CModelSPC : public Classifier
{
public:
  CModelSPC(const Mtx& objParams, int numStoredExemplars);

  virtual ~CModelSPC();

  virtual RequestResult handleRequest(fstring action,
                                      const Mtx& allModelParams,
                                      const MxWrapper& extraArgs);

private:
  virtual void computeDiffEv(const Mtx& objects,
                             Slice& modelParams, Mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

  const int itsNumStoredExemplars;
  const Mtx itsHiLo0;
  const Mtx itsHiLo1;
};

static const char vcid_cmodelspc_h[] = "$Header$";
#endif // !CMODELSPC_H_DEFINED
