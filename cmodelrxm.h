///////////////////////////////////////////////////////////////////////
//
// cmodelrxm.h
//
// Copyright (c) 2002-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
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
  CModelRxm(const mtx& objParams,
            TransferFunction transferFunc,
            int numStoredExemplars);

  virtual ~CModelRxm();

  virtual int numModelParams() const;

  virtual void loadModelParams(slice& modelParams);

  virtual const mtx& getStoredExemplars(Category cat);

protected:
  virtual int fillModelParamsBounds(mtx& bounds, int startRow) const;

private:
  mtx itsStored1;
  mtx itsStored2;

  const mtx itsHiLo1;
  const mtx itsHiLo2;
};

static const char vcid_cmodelrxm_h[] = "$Id$ $URL$";
#endif // !CMODELRXM_H_DEFINED
