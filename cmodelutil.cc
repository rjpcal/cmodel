///////////////////////////////////////////////////////////////////////
//
// cmodelutil.cc
//
// Copyright (c) 2002-2005
// Rob Peters <rjpeters at klab dot caltech dot edu>
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

void CModelUtil::clampRows(mtx& src, int firstrow, int nrows, const mtx& hilo)
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

mtx CModelUtil::getStoredExemplars(const slice& otherParams, int nstored,
                                   const mtx& hilo0, const mtx& hilo1)
{
DOTRACE("CModelUtil::getStoredExemplars");

  // reshape params into a stored exemplar matrix

  mtx storedExemplars =
    mtx(otherParams).as_shape(2*nstored, Classifier::DIM_OBJ_PARAMS);

  clampRows(storedExemplars, 0, nstored, hilo0);

  clampRows(storedExemplars, nstored, nstored, hilo1);

  return storedExemplars;
}

mtx CModelUtil::getHiLo(const mtx& src)
{
DOTRACE("CModelUtil::getHiLo");

  mtx result = mtx::zeros(2, src.ncols());

  for (int i = 0; i < src.ncols(); ++i)
    {
      result.at(0, i) = src.column(i).min();
      result.at(1, i) = src.column(i).max();
    }

  return result;
}

static const char vcid_cmodelutil_cc[] = "$Id$ $URL$";
#endif // !CMODELUTIL_CC_DEFINED
