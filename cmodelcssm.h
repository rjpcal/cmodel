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

#include "util/arrays.h"

class CModelCssm : public CModelExemplar {
private:
  const int itsNumTrainingExemplars;

  typedef const double* constDblPtr;

  fixed_block<constDblPtr> itsCat1;
  fixed_block<constDblPtr> itsCat2;

  // These will change depending on the most recent call to
  // findStoredExemplar
  double itsStored1[DIM_OBJ_PARAMS];
  double itsStored2[DIM_OBJ_PARAMS];

  const double* itsScaledWeights;

  // Scales the weights in place; weights is an input/output argument
  void scaleWeights(double* weights, int numRawWeights);

public:
  CModelCssm(const Rat& objParams,
				 const Rat& observedIncidence,
				 int numStoredExemplars);

  virtual ~CModelCssm();

  virtual void loadModelParams(Rat& modelParams);

  // The result of this function is only valid until the next call to
  // the function
  virtual const double* findStoredExemplar(Category cat, int n);
};

static const char vcid_cmodelcssm_h[] = "$Header$";
#endif // !CMODELCSSM_H_DEFINED
