///////////////////////////////////////////////////////////////////////
//
// cmodelutil.cc
//
// Copyright (c) 2002-2004 Rob Peters <rjpeters at klab dot caltech dot edu>
//
// created: Wed Jul 31 14:52:44 2002
// $Id$
//
///////////////////////////////////////////////////////////////////////

#ifndef CMODELUTIL_CC_DEFINED
#define CMODELUTIL_CC_DEFINED

#include "cmodel/cmodelutil.h"

#include "cmodel/classifier.h"

#include "mtx/mtx.h"

#include "util/trace.h"

void CModelUtil::clampRows(Mtx& src, int firstrow, int nrows, const Mtx& hilo)
{
DOTRACE("CModelUtil::clampRows");

  for (int c = 0; c < src.ncols(); ++c)
    {
      double clo = hilo.at(0, c);
      double chi = hilo.at(1, c);
      for (int r = 0; r < nrows; ++r)
        {
          double v = src.at(firstrow+r,c);
          if (v < clo) src.at(firstrow+r,c) = clo;
          if (v > chi) src.at(firstrow+r,c) = chi;
        }
    }
}

Mtx CModelUtil::getStoredExemplars(const Slice& otherParams, int nstored,
                                   const Mtx& hilo0, const Mtx& hilo1)
{
DOTRACE("CModelUtil::getStoredExemplars");

  // reshape params into a stored exemplar matrix

  Mtx storedExemplars =
    Mtx(otherParams).as_shape(2*nstored, Classifier::DIM_OBJ_PARAMS);

  clampRows(storedExemplars, 0, nstored, hilo0);

  clampRows(storedExemplars, nstored, nstored, hilo1);

  return storedExemplars;
}

Mtx CModelUtil::getHiLo(const Mtx& src)
{
DOTRACE("CModelUtil::getHiLo");

  Mtx result(2, src.ncols());

  for (int i = 0; i < src.ncols(); ++i)
    {
      result.at(0, i) = src.column(i).min();
      result.at(1, i) = src.column(i).max();
    }

  return result;
}

static const char vcid_cmodelutil_cc[] = "$Header$";
#endif // !CMODELUTIL_CC_DEFINED
