///////////////////////////////////////////////////////////////////////
//
// cmodelexemplar.cc
//
// Copyright (c) 1998-2000 Rob Peters rjpeters@klab.caltech.edu
//
// created: Fri Mar  9 14:32:31 2001
// written: Wed Mar 28 09:16:56 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELEXEMPLAR_CC_DEFINED
#define CMODELEXEMPLAR_CC_DEFINED

#include "cmodelexemplar.h"

#include "error.h"
#include "minivec.h"
#include "mtx.h"
#include "num.h"

#include "trace.h"

#include <cmath>

//
// This is just a functor that binds arguments to minkDist, so that
// we don't have to keep passing extra arguments to functions all the
// time in the time-critical inner loop.
//

class MinkowskiBinder {
public:
  MinkowskiBinder(MtxConstIter attWeights, MtxConstIter x2,
						 double r = 2.0, double r_inv = 0.5) :
	 itsAttWeights(attWeights),
	 itsX2(x2),
	 itsR(r),
	 itsRinv(r_inv)
  {}

  double minkDist(MtxConstIter x1) const
  {
	 double wt_sum = 0.0;
	 MtxConstIter wt = itsAttWeights;
	 MtxConstIter x2 = itsX2;

	 for (; wt.hasMore(); ++wt, ++x1, ++x2)
		{
		  wt_sum += (*wt) * pow( abs( *x1 - *x2), itsR);
		}
	 return pow(wt_sum, itsRinv);
  }

  // Specialized Minkowski distance for r==2
  double minkDist2(MtxConstIter x1) const;

private:
  const MtxConstIter itsAttWeights;
  const MtxConstIter itsX2;
  const double itsR;
  const double itsRinv;
};


///////////////////////////////////////////////////////////////////////
//
// CModelExemplar member definitions
//
///////////////////////////////////////////////////////////////////////

CModelExemplar::CModelExemplar(const Mtx& objParams,
										 const Mtx& observedIncidence,
										 int numStoredExemplars,
										 TransferFunction transferFunc) :
  Classifier(objParams, observedIncidence),
  itsNumTrainingExemplars(countCategory(objParams,0)),
  itsTraining1(itsNumTrainingExemplars, DIM_OBJ_PARAMS),
  itsTraining2(itsNumTrainingExemplars, DIM_OBJ_PARAMS),
  itsStored1Cache(numStoredExemplars, DIM_OBJ_PARAMS),
  itsStored2Cache(numStoredExemplars, DIM_OBJ_PARAMS),
  itsEvidence1Cache(numStoredExemplars, objParams.mrows()),
  itsEvidence2Cache(numStoredExemplars, objParams.mrows()),
  itsNumStoredExemplars(numStoredExemplars),
  itsTransferFunc(transferFunc)
{
  int num2 = countCategory(objParams, 1);

  if (itsNumTrainingExemplars != num2) {
	 throw ErrorWithMsg("the two categories must have the "
							  "same number of training exemplars");
  }

  // Find the category 1 and category 2 training exemplars
  int c1=0,c2=0;
  for (int i = 0; i < objParams.mrows(); ++i)
	 {
		if (int(objParams.at(i,0)) == 0)
		  {
			 itsTraining1.row(c1++) = exemplar(i);
		  }
		else if (int(objParams.at(i,0)) == 1)
		  {
			 itsTraining2.row(c2++) = exemplar(i);
		  }
	 }
}

CModelExemplar::~CModelExemplar()
{}

// Count the category training exemplars
int CModelExemplar::countCategory(const Mtx& params, int category) {
  int n = 0;
  for (int i = 0; i < params.mrows(); ++i)
	 {
		if (int(params.at(i,0)) == category)
		  ++n;
	 }
  return n;
}

//---------------------------------------------------------------------
//
// compute the minus loglikelihood for the constrained summed
// similarity model (cssm)
//
// based on the MATLAB function llcssm 
//
//     function ll = llcssm(modelParams, objParams, ...
//                          observedIncidence, numStoredExemplars)
//
//---------------------------------------------------------------------

void CModelExemplar::computeDiffEv(Slice& modelParams, Mtx& diffEvOut) {
DOTRACE("CModelExemplar::computeDiffEv");

  diffEvOut.setAll(0.0);

  //---------------------------------------------------------------------
  //
  // Set up the attentional weights.
  //

  Slice attWeights = modelParams.leftmost(DIM_OBJ_PARAMS);

  attWeights.apply(abs);

  Slice otherParams = modelParams.rightmost(modelParams.nelems()-
														  (DIM_OBJ_PARAMS+2));

  loadModelParams(otherParams);

  //---------------------------------------------------------------------
  //
  // Compute diffEvidence, a matrix of differences of summed similarities.
  //

  const double minkPower = 2.0;
  const double minkPowerInv = 1.0/minkPower;

  const Mtx& stored1 = getStoredExemplars(CAT1);
  const Mtx& stored2 = getStoredExemplars(CAT2);

  MtxConstIter attWts = attWeights.begin();

  minivec<MtxConstIter> exemplars;

  for (int yy = 0; yy < numAllExemplars(); ++yy) {
	 exemplars.push_back(exemplar(yy).begin());
  }

  const MtxIter diffEv = diffEvOut.colIter(0);

  for (int x = 0; x < itsNumStoredExemplars; ++x) {
	 DOTRACE("minkowski loop");
	 MinkowskiBinder binder1(attWts, stored1.rowIter(x),
									 minkPower, minkPowerInv);
	 MinkowskiBinder binder2(attWts, stored2.rowIter(x),
									 minkPower, minkPowerInv);

	 const MtxIter distrust1 = itsEvidence1Cache.rowIter(x);
	 const MtxIter distrust2 = itsEvidence2Cache.rowIter(x);

  	 bool compute1 = true;
  	 bool compute2 = true;

//  	 bool compute1 = (itsStored1Cache.row(x) != stored1.row(x));
//  	 bool compute2 = (itsStored2Cache.row(x) != stored2.row(x));

	 if (compute1) {
		DOTRACE("compute1");
		itsStored1Cache.row(x) = stored1.row(x);

		int y = 0;
		for (MtxIter iter1 = distrust1;
			  iter1.hasMore();
			  ++y, ++iter1) {
		  if (minkPower == 2.0) {
			 if (EXP_DECAY == itsTransferFunc) {
				*iter1 = Num::fastexp7(-binder1.minkDist2(exemplars[y]));
			 }
			 else if (LINEAR_DECAY == itsTransferFunc) {
				*iter1 = -binder1.minkDist2(exemplars[y]);
			 }
		  }
		  else {
			 if (EXP_DECAY == itsTransferFunc) {
				*iter1 = Num::fastexp7(-binder1.minkDist(exemplars[y]));
			 }
			 else if (LINEAR_DECAY == itsTransferFunc) {
				*iter1 = -binder1.minkDist(exemplars[y]);
			 }
		  }
		}
	 }

	 if (compute2) {
		DOTRACE("compute1");
		itsStored2Cache.row(x) = stored2.row(x);

		int y = 0;
		for (MtxIter iter2 = distrust2;
			  iter2.hasMore();
			  ++y, ++iter2) {
		  if (minkPower == 2.0) {
			 if (EXP_DECAY == itsTransferFunc) {
				*iter2 = Num::fastexp7(-binder2.minkDist2(exemplars[y]));
			 }
			 else if (LINEAR_DECAY == itsTransferFunc) {
				*iter2 = -binder2.minkDist2(exemplars[y]);
			 }
		  }
		  else {
			 if (EXP_DECAY == itsTransferFunc) {
				*iter2 = Num::fastexp7(-binder2.minkDist(exemplars[y]));
			 }
			 else if (LINEAR_DECAY == itsTransferFunc) {
				*iter2 = -binder2.minkDist(exemplars[y]);
			 }
		  }
		}
	 }

  }

  for (int x = 0; x < itsNumStoredExemplars; ++x) {

	 const MtxIter distrust1 = itsEvidence1Cache.rowIter(x);
	 const MtxIter distrust2 = itsEvidence2Cache.rowIter(x);

	 for (MtxIter iter1 = distrust1, iter2 = distrust2, diff = diffEv;
			iter1.hasMore();
			++iter1, ++iter2, ++diff) {
		*diff += *iter1 - *iter2;
	 }
  }
}

double CModelExemplar::computeSigmaNoise(double rawSigma) const
{
  return rawSigma * sqrt(itsNumStoredExemplars*2.0);
}

void CModelExemplar::loadModelParams(Slice& modelParams) {}

#define USE_ASM

double MinkowskiBinder::minkDist2(MtxConstIter x1) const
{
#ifndef USE_ASM

  double wt_sum = 0.0;
  MtxConstIter wt = itsAttWeights;
  MtxConstIter x2 = itsX2;

  for (; wt.hasMore(); ++wt, ++x1, ++x2)
    {
      wt_sum += (*wt) * ((*x1) - (*x2)) * ((*x1) - (*x2));
    }
  return sqrt(wt_sum);

#else

// layout of *this
//      12 bytes -- wts
// +0        4 bytes -- data pointer
// +4        4 bytes -- stride
// +8        4 bytes -- stop pointer
//      12 bytes -- x2
// +12       4 bytes -- data pointer
// +16       4 bytes -- stride
// +20       4 bytes -- stop pointer
// +24  8 bytes  -- itsR
// +32  8 bytes  -- itsRinv

asm(
    // "pushl %ebp" assumed
    // "movl %esp,%ebp" assumed

    "subl $4,%esp;"             // {1c} 4b on stack
    "movl 8(%ebp),%edx;"        // {1c} copy the "this" pointer to edx

    "fldz;"                     // {2c} initialize the result to 0.0

    "pushl %edi;"					  // {1c}
    "pushl %esi;"					  // {1c}
    "pushl %ebx;"					  // {1c}

#define WTS_REG "%eax"
    "movl (%edx),"WTS_REG";"    // {1c} copy wts.data to WTS_REG

    "movl 8(%edx),%ebx;"        // {1c} copy wts.stop to ebx

#define X2_REG "%ecx"
    "movl 12(%edx),"X2_REG";"   // {1c} copy x2.data to X2_REG

    "cmpl %ebx,"WTS_REG";"		  // {1c}
    "jae cleanup;"              // {3/1c} jump if (wts.data >= wts.stop)

    "movl 4(%edx),%edi;"        // {1c} copy wts.stride to edi
    "sall $3,%edi;"             // {2c} wts.stride *= 8 (which is sizeof(double))
    "movl %edi,-4(%ebp);"       // {1c} copy wts.stride to stack-4

    "movl 16(%ebp),%edi;"       // {1c} copy x1.stride to edi
    "sall $3,%edi;"             // {2c} x1.stride *= 8

    "movl 16(%edx),%esi;"       // {1c} copy x2.stride to esi

    "sall $3,%esi;"             // {2c} x2.stride *= 8
    ".p2align 4,,7;"

#define X1_REG "%edx"

    "movl 12(%ebp),"X1_REG";"   // {1c} copy x1.data to edx

"sum_loop:;"
    "fldl ("X1_REG");"          // {2c} top = *x1.data
    "fsubl ("X2_REG");"         // {3/1c} top -= *x2.data
    "fldl ("WTS_REG");"         // {2c} top = *wts.data
    "fmul %st(1),%st;"          // {3/1c} top *= (*x1.data - *x2.data)

	 "addl -4(%ebp),"WTS_REG";"  // {2c} wts.data += wts.stride
  	 "addl %edi,"X1_REG";"       // {1c} x1.data += x1.stride
  	 "addl %esi,"X2_REG";"       // {1c} x2.data += x2.stride

    "fmulp %st,%st(1);"         // {3/1c} top *= (*x1.data - *x2.data) and pop st(0)
    "faddp %st,%st(1);"         // {3/1c} add term to wt_sum and pop st(0)

    "cmpl %ebx,"WTS_REG";"		  // {1c}
    "jb sum_loop;"              // {3/1c} loop if (wts.data < wts.stop)

	 "fsqrt;"                    // {70c} st(0) = sqrt(wt_sum)

"cleanup:;"
    "leal -16(%ebp),%esp;"		  // {1c}
    "popl %ebx;"					  // {4c}
    "popl %esi;"					  // {4c}
    "popl %edi;"					  // {4c}
    // "movl %ebp,%esp;" assumed
    // "popl %ebp;" assumed
    // "ret;" assumed
    ); // end asm
#endif
}

static const char vcid_cmodelexemplar_cc[] = "$Header$";
#endif // !CMODELEXEMPLAR_CC_DEFINED
