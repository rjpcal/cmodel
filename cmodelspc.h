///////////////////////////////////////////////////////////////////////
//
// cmodelspc.h
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Mon Feb  4 13:59:00 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELSPC_H_DEFINED
#define CMODELSPC_H_DEFINED

#include "cmodel/classifier.h"

class CModelSPC : public Classifier
{
public:
  CModelSPC(const Mtx& objParams, int numStoredExemplars);

  virtual ~CModelSPC();

  virtual int numModelParams() const;

  virtual RequestResult handleRequest(fstring action,
                                      const Mtx& allModelParams,
                                      const MxWrapper& extraArgs);

protected:
  virtual int fillModelParamsBounds(Mtx& bounds, int startRow) const;

private:
  virtual void computeDiffEv(const Mtx& objects,
                             slice& modelParams, Mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

  const int itsNumStoredExemplars;
  const Mtx itsHiLo0;
  const Mtx itsHiLo1;
};

static const char vcid_cmodelspc_h[] = "$Header$";
#endif // !CMODELSPC_H_DEFINED
