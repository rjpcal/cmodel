///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:40:21 2001
// written: Thu Apr 26 18:41:27 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_H_DEFINED
#define CMODELGCM_H_DEFINED

#include "cmodelexemplar.h"

class CModelGcm : public CModelExemplar {
private:
  virtual const Mtx& getStoredExemplars(Category cat);

public:
  CModelGcm(const Mtx& objParams,
				TransferFunction transferFunc = EXP_DECAY);

  virtual ~CModelGcm();
};

static const char vcid_cmodelgcm_h[] = "$Header$";
#endif // !CMODELGCM_H_DEFINED
