///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:24:41 2001
// written: Thu Apr 26 18:40:56 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_H_DEFINED
#define CMODELCSSM_H_DEFINED

#include "cmodelexemplar.h"

class CModelCssm : public CModelExemplar {
private:
  Mtx itsStored1;
  Mtx itsStored2;

  Mtx itsCachedRawWts1;
  Mtx itsCachedRawWts2;

public:
  CModelCssm(const Mtx& objParams,
				 TransferFunction transferFunc,
				 int numStoredExemplars);

  virtual ~CModelCssm();

  virtual void loadModelParams(Slice& modelParams);

  virtual const Mtx& getStoredExemplars(Category cat);
};

static const char vcid_cmodelcssm_h[] = "$Header$";
#endif // !CMODELCSSM_H_DEFINED
