///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.h
//
// Copyright (c) 2001-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Fri Mar  9 14:40:21 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_H_DEFINED
#define CMODELGCM_H_DEFINED

#include "cmodel/cmodelexemplar.h"

class CModelGcm : public CModelExemplar
{
private:
  virtual const mtx& getStoredExemplars(Category cat);

public:
  CModelGcm(const mtx& objParams,
            TransferFunction transferFunc = EXP_DECAY);

  virtual ~CModelGcm();
};

static const char vcid_cmodelgcm_h[] = "$Id$ $URL$";
#endif // !CMODELGCM_H_DEFINED
