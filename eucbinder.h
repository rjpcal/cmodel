///////////////////////////////////////////////////////////////////////
//
// eucbinder.h
//
// Copyright (c) 2002-2002 Rob Peters rjpeters@klab.caltech.edu
//
// created: Thu Feb  7 16:18:49 2002
// written: Wed Feb 20 18:10:31 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef EUCBINDER_H_DEFINED
#define EUCBINDER_H_DEFINED

#if defined(NO_EXTERNAL_INCLUDE_GUARDS) || !defined(MTX_H_DEFINED)
#include "mtx/mtx.h"
#endif

//
// This is just a functor that binds arguments to eucDist, so that we don't
// have to keep passing extra arguments to functions all the time in the
// time-critical inner loop.
//
// In addition, this class is templatized on the number of elements in each
// vector, so that loops can be better optimized at compile-time. Therefore
// note that no range-checking will be done at run-time.
//

template <unsigned int Dim>
class EuclideanBinder
{
public:
  EuclideanBinder(MtxConstIter attWeights, MtxConstIter x2) :
    itsAttWeightsStop(&itsAttWeights[0] + Dim)
  {
    double* wtptr = &itsAttWeights[0];
    double* x2ptr = &itsX2[0];

    for (; wtptr < itsAttWeightsStop; ++wtptr, ++x2ptr, ++attWeights, ++x2)
      {
        *wtptr = *attWeights;
        *x2ptr = *x2;
      }
  }

  double eucDist(MtxConstIter x1) const
  {
    double wt_sum = 0.0;

    const double* wt = &itsAttWeights[0];
    const double* x2 = &itsX2[0];

    for (; wt < itsAttWeightsStop; ++wt, ++x2, ++x1)
      {
        wt_sum += (*wt) * square((*x1) - (*x2));
      }
    return sqrt(wt_sum);
  }

private:
  static double square(double x) { return x * x; }

  double itsAttWeights[Dim];
  double* const itsAttWeightsStop;
  double itsX2[Dim];
};

static const char vcid_eucbinder_h[] = "$Header$";
#endif // !EUCBINDER_H_DEFINED
