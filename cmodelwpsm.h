///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.h
//
// Copyright (c) 2001-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 17:35:37 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELWPSM_H_DEFINED
#define CMODELWPSM_H_DEFINED

#include "cmodel/cmodelexemplar.h"

class CModelWpsm : public CModelExemplar
{
private:
  mtx itsPrototype1;
  mtx itsPrototype2;

  virtual const mtx& getStoredExemplars(Category cat);

public:
  CModelWpsm(const mtx& objParams,
             TransferFunction transferFunc = EXP_DECAY);

  virtual ~CModelWpsm();
};

static const char vcid_cmodelwpsm_h[] = "$Id$ $URL$";
#endif // !CMODELWPSM_H_DEFINED
