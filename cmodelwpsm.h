///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.h
//
// Copyright (c) 1998-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 17:35:37 2001
// written: Thu Feb  7 13:50:56 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELWPSM_H_DEFINED
#define CMODELWPSM_H_DEFINED

#include "cmodelexemplar.h"

class CModelWpsm : public CModelExemplar
{
private:
  Mtx itsPrototype1;
  Mtx itsPrototype2;

  virtual const Mtx& getStoredExemplars(Category cat);

public:
  CModelWpsm(const Mtx& objParams,
             TransferFunction transferFunc = EXP_DECAY);

  virtual ~CModelWpsm();
};

static const char vcid_cmodelwpsm_h[] = "$Header$";
#endif // !CMODELWPSM_H_DEFINED
