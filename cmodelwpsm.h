///////////////////////////////////////////////////////////////////////
//
// cmodelwpsm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 17:35:37 2001
// written: Tue Mar 13 18:02:28 2001
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
