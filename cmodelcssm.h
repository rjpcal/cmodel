///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.h
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Thu Mar  8 16:24:41 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_H_DEFINED
#define CMODELCSSM_H_DEFINED

#include "cmodel/cmodelexemplar.h"

class CModelCssm : public CModelExemplar
{
public:
  CModelCssm(const mtx& objParams,
             TransferFunction transferFunc,
             int numStoredExemplars);

  virtual ~CModelCssm();

  virtual int numModelParams() const;

  virtual void loadModelParams(slice& modelParams);

  virtual const mtx& getStoredExemplars(Category cat);

protected:
  virtual int fillModelParamsBounds(mtx& bounds, int startRow) const;

private:
  mtx itsStored1;
  mtx itsStored2;

  mtx itsCachedRawWts1;
  mtx itsCachedRawWts2;
};

static const char vcid_cmodelcssm_h[] = "$Header$";
#endif // !CMODELCSSM_H_DEFINED
