///////////////////////////////////////////////////////////////////////
//
// minkbinder.h
//
// Copyright (c) 2001-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Mon Jul  9 13:59:12 2001
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef MINKBINDER_H_DEFINED
#define MINKBINDER_H_DEFINED

#include "mtx/mtx.h"

#include <cmath>

//
// This is just a functor that binds arguments to minkDist, so that
// we don't have to keep passing extra arguments to functions all the
// time in the time-critical inner loop.
//

class MinkowskiBinder
{
public:
  MinkowskiBinder(mtx_const_iter attWeights, mtx_const_iter x2,
                  double r = 2.0, double r_inv = 0.5) :
    itsAttWeights(attWeights),
    itsX2(x2),
    itsR(r),
    itsRinv(r_inv)
  {}

  double minkDist(mtx_const_iter x1) const
  {
    double wt_sum = 0.0;
    mtx_const_iter wt = itsAttWeights;
    mtx_const_iter x2 = itsX2;

    for (; wt.has_more(); ++wt, ++x1, ++x2)
      {
        wt_sum += (*wt) * pow( std::abs( *x1 - *x2), itsR);
      }
    return pow(wt_sum, itsRinv);
  }

  // Specialized Minkowski distance for r==2
  double minkDist2(mtx_const_iter x1) const;

private:
  const mtx_const_iter itsAttWeights;
  const mtx_const_iter itsX2;
  const double itsR;
  const double itsRinv;
};

static const char vcid_minkbinder_h[] = "$Header$";
#endif // !MINKBINDER_H_DEFINED
