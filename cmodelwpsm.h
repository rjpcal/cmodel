///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 17:35:37 2001
// written: Mon Mar 12 16:44:59 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELWPSM_H_DEFINED
#define CMODELWPSM_H_DEFINED

#include "cmodelexemplar.h"

class CModelWpsm : public CModelExemplar {
private:
  double itsPrototype1[DIM_OBJ_PARAMS];
  double itsPrototype2[DIM_OBJ_PARAMS];

  virtual void loadModelParams(Mtx& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual ConstSlice findStoredExemplar(Category cat, int n);

public:
  CModelWpsm(const Mtx& objParams,
				 const Mtx& observedIncidence);

  virtual ~CModelWpsm();
};

static const char vcid_cmodelwpsm_h[] = "$Header$";
#endif // !CMODELWPSM_H_DEFINED
