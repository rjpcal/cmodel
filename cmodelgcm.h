///////////////////////////////////////////////////////////////////////
//
// cmodelgcm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:40:21 2001
// written: Fri Mar  9 16:51:17 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELGCM_H_DEFINED
#define CMODELGCM_H_DEFINED

#include "cmodelexemplar.h"

#include "util/arrays.h"

class CModelGcm : public CModelExemplar {
private:
  const int itsNumTrainingExemplars;

  typedef const double* constDblPtr;

  fixed_block<constDblPtr> itsCat1;
  fixed_block<constDblPtr> itsCat2;

  virtual void loadModelParams(Rat& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual Slice findStoredExemplar(Category cat, int n);

public:
  CModelGcm(const Rat& objParams,
				const Rat& observedIncidence);

  virtual ~CModelGcm();
};

static const char vcid_cmodelgcm_h[] = "$Header$";
#endif // !CMODELGCM_H_DEFINED
