///////////////////////////////////////////////////////////////////////
//
// cmodelcssm.h
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Mar  8 16:24:41 2001
// written: Fri Mar  9 14:31:51 2001
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
  double itsStored1[DIM_OBJ_PARAMS];
  double itsStored2[DIM_OBJ_PARAMS];

//    const double* itsScaledWeights;

  ConstSlice itsScaledWeights;

  // Scales the weights in place; weights is an input/output argument
  void scaleWeights(Slice& weights);

public:
  CModelCssm(const Mtx& objParams,
				 const Mtx& observedIncidence,
				 int numStoredExemplars);

  virtual ~CModelCssm();

  virtual void loadModelParams(Mtx& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual ConstSlice findStoredExemplar(Category cat, int n);
};

static const char vcid_cmodelcssm_h[] = "$Header$";
#endif // !CMODELCSSM_H_DEFINED
