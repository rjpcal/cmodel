///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:40:21 2001
// written: Mon Mar 12 16:44:54 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_H_DEFINED
#define CMODELGCM_H_DEFINED

#include "cmodelexemplar.h"

class CModelGcm : public CModelExemplar {
private:
  virtual void loadModelParams(Mtx& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual ConstSlice findStoredExemplar(Category cat, int n);

public:
  CModelGcm(const Mtx& objParams,
				const Mtx& observedIncidence);

  virtual ~CModelGcm();
};

static const char vcid_cmodelgcm_h[] = "$Header$";
#endif // !CMODELGCM_H_DEFINED
