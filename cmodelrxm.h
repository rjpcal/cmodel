///////////////////////////////////////////////////////////////////////
//
// cmodelrxm.h
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Wed Jul 31 14:47:13 2002
// written: Wed Jul 31 14:47:13 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELRXM_H_DEFINED
#define CMODELRXM_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(CMODELEXEMPLAR_H_DEFINED)
#include "cmodel/cmodelexemplar.h"
#endif

class CModelRxm : public CModelExemplar
{
public:
  CModelRxm(const Mtx& objParams,
	    TransferFunction transferFunc,
	    int numStoredExemplars);

  virtual ~CModelRxm();

  virtual void loadModelParams(Slice& modelParams);

  virtual const Mtx& getStoredExemplars(Category cat);

private:
  Mtx itsStored1;
  Mtx itsStored2;

  const Mtx itsHiLo1;
  const Mtx itsHiLo2;
};

static const char vcid_cmodelrxm_h[] = "$Header$";
#endif // !CMODELRXM_H_DEFINED
