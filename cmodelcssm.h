///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:24:41 2001
// written: Tue Sep 28 14:26:46 2004
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_H_DEFINED
#define CMODELCSSM_H_DEFINED

#include "cmodel/cmodelexemplar.h"

class CModelCssm : public CModelExemplar
{
public:
  CModelCssm(const Mtx& objParams,
             TransferFunction transferFunc,
             int numStoredExemplars);

  virtual ~CModelCssm();

  virtual int numModelParams() const;

  virtual void loadModelParams(Slice& modelParams);

  virtual const Mtx& getStoredExemplars(Category cat);

protected:
  virtual int fillModelParamsBounds(Mtx& bounds, int startRow) const;

private:
  Mtx itsStored1;
  Mtx itsStored2;

  Mtx itsCachedRawWts1;
  Mtx itsCachedRawWts2;
};

static const char vcid_cmodelcssm_h[] = "$Header$";
#endif // !CMODELCSSM_H_DEFINED
