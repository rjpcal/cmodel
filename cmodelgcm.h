///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:40:21 2001
// written: Tue Sep 28 14:26:54 2004
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_H_DEFINED
#define CMODELGCM_H_DEFINED

#include "cmodel/cmodelexemplar.h"

class CModelGcm : public CModelExemplar
{
private:
  virtual const Mtx& getStoredExemplars(Category cat);

public:
  CModelGcm(const Mtx& objParams,
            TransferFunction transferFunc = EXP_DECAY);

  virtual ~CModelGcm();
};

static const char vcid_cmodelgcm_h[] = "$Header$";
#endif // !CMODELGCM_H_DEFINED
