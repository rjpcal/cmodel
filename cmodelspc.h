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
  CModelSPC(const mtx& objParams, int numStoredExemplars);

  virtual ~CModelSPC();

  virtual int numModelParams() const;

  virtual RequestResult handleRequest(fstring action,
                                      const mtx& allModelParams,
                                      const MxWrapper& extraArgs);

protected:
  virtual int fillModelParamsBounds(mtx& bounds, int startRow) const;

private:
  virtual void computeDiffEv(const mtx& objects,
                             slice& modelParams, mtx& diffEvOut);

  virtual double computeSigmaNoise(double rawSigma) const;

  const int itsNumStoredExemplars;
  const mtx itsHiLo0;
  const mtx itsHiLo1;
};

static const char vcid_cmodelspc_h[] = "$Header$";
#endif // !CMODELSPC_H_DEFINED
