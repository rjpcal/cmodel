///////////////////////////////////////////////////////////////////////
//
// cmodelrxm.h
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Wed Jul 31 14:47:13 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELRXM_H_DEFINED
#define CMODELRXM_H_DEFINED

#include "cmodel/cmodelexemplar.h"

class CModelRxm : public CModelExemplar
{
public:
  CModelRxm(const Mtx& objParams,
            TransferFunction transferFunc,
            int numStoredExemplars);

  virtual ~CModelRxm();

  virtual int numModelParams() const;

  virtual void loadModelParams(Slice& modelParams);

  virtual const Mtx& getStoredExemplars(Category cat);

protected:
  virtual int fillModelParamsBounds(Mtx& bounds, int startRow) const;

private:
  Mtx itsStored1;
  Mtx itsStored2;

  const Mtx itsHiLo1;
  const Mtx itsHiLo2;
};

static const char vcid_cmodelrxm_h[] = "$Header$";
#endif // !CMODELRXM_H_DEFINED
