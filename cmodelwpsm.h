///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.h
//
// Copyright (c) 1998-2001 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 17:35:37 2001
// written: Fri Apr  6 10:27:20 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELWPSM_H_DEFINED
#define CMODELWPSM_H_DEFINED

#include "cmodelexemplar.h"

class CModelWpsm : public CModelExemplar {
private:
  Mtx itsPrototype1;
  Mtx itsPrototype2;

  virtual const Mtx& getStoredExemplars(Category cat);

public:
  CModelWpsm(const Mtx& objParams,
				 const Mtx& observedIncidence,
				 TransferFunction transferFunc = EXP_DECAY);

  virtual ~CModelWpsm();
};

static const char vcid_cmodelwpsm_h[] = "$Header$";
#endif // !CMODELWPSM_H_DEFINED
