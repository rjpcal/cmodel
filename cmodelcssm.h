///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:24:41 2001
// written: Tue Mar 13 18:01:31 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELCSSM_H_DEFINED
#define CMODELCSSM_H_DEFINED

#include "cmodelexemplar.h"

class CModelCssm : public CModelExemplar {
private:
  // These will change depending on the most recent call to
  // findStoredExemplar
  Mtx itsStored1;
  Mtx itsStored2;

  Mtx itsScaledWeights;

public:
  CModelCssm(const Mtx& objParams,
				 const Mtx& observedIncidence,
				 int numStoredExemplars);

  virtual ~CModelCssm();

  virtual void loadModelParams(Slice& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual ConstSlice findStoredExemplar(Category cat, int n);
};

static const char vcid_cmodelcssm_h[] = "$Header$";
#endif // !CMODELCSSM_H_DEFINED
