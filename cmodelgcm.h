///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.h
//
// Copyright (c) 2001-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:40:21 2001
// written: Thu Feb 14 11:54:37 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_H_DEFINED
#define CMODELGCM_H_DEFINED

#include "cmodelexemplar.h"

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
